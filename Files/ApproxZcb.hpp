//
//  ApproxZcb.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 10/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef ApproxZcb_hpp
#define ApproxZcb_hpp

#include <stdio.h>
#include <ql/models/shortrate/twofactormodel.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>


inline void approxZcb() {
    Date anchor_date(26,May,2009);
    Settings::instance().evaluationDate() = anchor_date;
    
    Actual360 zero_day_count;
    TARGET calendar;
    
    std::vector<Date> the_dates;
    std::vector<Rate> the_rates;
    
    std::string file_name = "/Users/mmencke/Speciale_Data/table12-1.csv";
    
    
    Matrix data = readCsv(file_name);
    /*
     std::cout << std::endl;
     printMatrix(data);
     std::cout << std::endl;
     */
    for(int i=0;i<data.rows();i++) {
        the_dates.push_back(Date(data[i][0]));
        the_rates.push_back(data[i][1]);
    }
    
    
    ext::shared_ptr<YieldTermStructure> the_curve_ptr(new InterpolatedZeroCurve<Cubic>(the_dates,the_rates, zero_day_count,calendar));
    
    RelinkableHandle<YieldTermStructure> the_handle;
    
    the_handle.linkTo(the_curve_ptr);
    
    ext::shared_ptr<G2> g2_model(new G2(the_handle));
    Array g2_params={0.0183149,0.0106535,0.0183039,0.0110041,-0.729877};
    g2_model->setParams(g2_params);
    
    
    ext::shared_ptr<TwoFactorModel::ShortRateDynamics> g2_dynamics =g2_model->dynamics();
    
    ext::shared_ptr<StochasticProcess> g2_process = g2_dynamics->process();
    

    



    
    unsigned long seed =1;
    
    MersenneTwisterUniformRng uniformGenerator(seed);

    InverseCumulativeRng<MersenneTwisterUniformRng,InverseCumulativeNormal> rng(uniformGenerator);
    
    

    Natural n = 10000;
    Real dt=1.0/52.0;
    Real sum_approx=0,sum_true=0;
    
    
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff_approx= end-start;
    auto diff_true=end-start;
    
    
    for(int k=0;k<n;k++) {
        Time t=0.0;
        Real zcb_approx=1;
        Real zcb_true=1;
        
        Array x=g2_process->initialValues();

        
        for(int i=0;i<520;i++) {
            Array dw(x.size());
            for(int j=0;j<x.size();j++) {
                dw[j]=rng.next().value;
            }
            x=g2_process->evolve(t,x,dt,dw);
            
            start = std::chrono::steady_clock::now();
            zcb_approx*=std::exp(-g2_dynamics->shortRate(t,x[0],x[1])*dt);
            end = std::chrono::steady_clock::now();
            diff_approx+= end-start;
            
            start = std::chrono::steady_clock::now();
            zcb_true*=g2_model->discountBond(t,t+dt,x);
            end = std::chrono::steady_clock::now();
            diff_true+=end-start;
            t+=dt;
        }
            sum_approx+=zcb_approx;
            sum_true+=zcb_true;
    }
    sum_approx/=n;
    sum_true/=n;
    
    std::cout << sum_approx << ", " << sum_true << std::endl;
    std::cout << "Approx = " << std::chrono::duration <double, std::milli>(diff_approx).count() << " ms" << std::endl;
    std::cout << "True = " << std::chrono::duration <double, std::milli>(diff_true).count() << " ms" << std::endl;
    
    
    Matrix rho(3,3);
    rho[0][0]=rho[1][1]=rho[2][2]=1;
    rho[0][1]=rho[1][0]=0;
    rho[0][2]=rho[2][0]=0.5;
    rho[1][2]=rho[2][1]=-0.5;
    
    std::cout << std::endl;
    printMatrix(rho);
    std::cout << std::endl;
    
    Matrix rho_chol=CholeskyDecomposition(rho);
    
    std::cout << std::endl;
    printMatrix(rho_chol);
    std::cout << std::endl;
    
    Matrix correlated_normals(n,3);
    
    
    
    for(int k=0;k<n;k++) {
        Array dw(3);
        for(int j=0;j<3;j++) {
            dw[j]=rng.next().value;
        }
        Array dw_cor(3);
        dw_cor=rho_chol*dw;
        
        for(int j=0;j<3;j++) {
            correlated_normals[k][j]=dw_cor[j];
        }
    }
    
    writeCsv(correlated_normals,"correlated_normals.csv");
}

#endif /* ApproxZcb_hpp */
