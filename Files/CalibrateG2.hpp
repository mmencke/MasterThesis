//
//  CalibrateG2.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 24/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef CalibrateG2_hpp
#define CalibrateG2_hpp

#include "utility.hpp"

#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>

#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>

#include <ql/pricingengines/swaption/g2swaptionengine.hpp>
//#include <ql/pricingengines/swaption/treeswaptionengine.hpp>
#include <ql/models/shortrate/twofactormodels/g2.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/endcriteria.hpp>

#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/pricingengines/swaption/jamshidianswaptionengine.hpp>

#include <boost/format.hpp>

using namespace QuantLib;



void calibrateG2() {
    
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
    
    
    Matrix term_structure(52*30,2);
    
    double t=0;
    double dt=(1/(double)52);
    for(int i=0;i<52*30;i++) {
        t=t+dt;
        
        term_structure[i][0]=t;
        term_structure[i][1]=the_curve_ptr->zeroRate(t, Continuous);
    }
    
    
    writeCsv(term_structure, "term_structure.csv");

    std::vector<Period> tenors;
    std::vector<Period> expiries;
    
    for(int i=1;i<=10;i++) tenors.push_back(i*Years);
    //for(int i=1;i<=10;i++) std::cout << tenors[i-1] << std::endl;
    
    for(int i=1;i<=5;i++) expiries.push_back(i*Years);
    expiries.push_back(7*Years);
    expiries.push_back(10*Years);
    
    
    Matrix vol_surface = readCsv("/Users/mmencke/Speciale_Data/table12-2.csv");
    
    /*
    std::cout << std::endl;
    printMatrix(vol_surface);
    std::cout << std::endl;
    */
    ext::shared_ptr<IborIndex> euribor_index(new Euribor3M(the_handle));
    
    Thirty360 fixed_dc;
    Period fixed_tenor=1*Years;
    Actual360 float_dc;
    
    std::vector<ext::shared_ptr<SwaptionHelper>> swaption_helper_vec;
    
    ext::shared_ptr<G2> g2_ptr(new G2(the_handle));
    //ext::shared_ptr<PricingEngine> the_engine_ptr(new TreeSwaptionEngine(g2_ptr,25));
    ext::shared_ptr<PricingEngine> the_engine_ptr(new G2SwaptionEngine(g2_ptr,4, 100));//±4 standard deviations wide with 100 intervals
    
    for(int i=0;i<expiries.size();i++) {
        for(int j=0;j<tenors.size();j++) {
            ext::shared_ptr<SwaptionHelper> temp_helper_ptr(new SwaptionHelper(expiries[i],tenors[j],Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(vol_surface[i][j]))), euribor_index, fixed_tenor,fixed_dc,float_dc,  the_handle));
            
            temp_helper_ptr->setPricingEngine(the_engine_ptr);
            //std::cout << "Model value=" <<temp_helper_ptr->modelValue() << std::endl;
            swaption_helper_vec.push_back(temp_helper_ptr);
        }
    }
    
 
    LevenbergMarquardt  optimisation_method(1.0e-8,1.0e-8,1.0e-8);//epsfcn,xtol,gtol
    EndCriteria end_criteria(1000, 100, 1e-6, 1e-8, 1e-8);//maxIterations,  maxStationaryStateIterations,  rootEpsilon,  functionEpsilon,  gradientNormEpsilon
   
    
    std::vector<ext::shared_ptr<CalibrationHelper>> temp_helper_vec;
    for(int i=0; i<swaption_helper_vec.size();i++) {
        temp_helper_vec.push_back(swaption_helper_vec[i]);
    }
    g2_ptr->calibrate(temp_helper_vec, optimisation_method, end_criteria);

    std::cout << "G2 Params" << std::endl;
    std::cout <<g2_ptr->params()[0]<<std::endl;
    std::cout <<g2_ptr->params()[1]<<std::endl;
    std::cout <<g2_ptr->params()[2]<<std::endl;
    std::cout <<g2_ptr->params()[3]<<std::endl;
    std::cout <<g2_ptr->params()[4]<<std::endl;
    std::cout << std::endl;

    Array g2_params={0.0183149,0.0106535,0.0183039,0.0110041,-0.729877};
    g2_ptr->setParams(g2_params);
    
    Matrix calibrated_g2(70,5);
    
    int n=expiries.size();
    int m=tenors.size();
    
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            calibrated_g2[i*m+j][0]=swaption_helper_vec[i*m+j]->marketValue();
            calibrated_g2[i*m+j][1]=swaption_helper_vec[i*m+j]->modelValue();
        }
    }
    
    Array g2_brigo_params = {0.0002,0.0080,7.6630,0.0182,-0.9734};
    
    g2_ptr->setParams(g2_brigo_params);

    
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            calibrated_g2[i*m+j][2]=swaption_helper_vec[i*m+j]->modelValue();
        }
    }
    
    
    ext::shared_ptr<HullWhite> hull_white(new HullWhite(the_handle));
    //ext::shared_ptr<PricingEngine> the_engine_ptr(new TreeSwaptionEngine(g2_ptr,25));
    ext::shared_ptr<PricingEngine> hw_engine(new JamshidianSwaptionEngine(hull_white,the_handle));
    
    
    temp_helper_vec.clear();
    //std::vector<ext::shared_ptr<CalibrationHelper>> temp_helper_vec;
    for(int i=0; i<swaption_helper_vec.size();i++) {
        swaption_helper_vec[i]->setPricingEngine(hw_engine);
        temp_helper_vec.push_back(swaption_helper_vec[i]);
    }
    hull_white->calibrate(temp_helper_vec, optimisation_method, end_criteria);
    
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            calibrated_g2[i*m+j][3]=swaption_helper_vec[i*m+j]->modelValue();
        }
    }
    std::cout << "Hull-White Params" << std::endl;
    std::cout <<hull_white->params()[0]<<std::endl;
    std::cout <<hull_white->params()[1]<<std::endl;
    std::cout << std::endl;
    
    
    ext::shared_ptr<Vasicek> vasicek(new Vasicek);
    ext::shared_ptr<PricingEngine> vasicek_engine(new JamshidianSwaptionEngine(vasicek,the_handle));
    
    temp_helper_vec.clear();
    //std::vector<ext::shared_ptr<CalibrationHelper>> temp_helper_vec;
    for(int i=0; i<swaption_helper_vec.size();i++) {
        swaption_helper_vec[i]->setPricingEngine(vasicek_engine);
        temp_helper_vec.push_back(swaption_helper_vec[i]);
    }
    vasicek->calibrate(temp_helper_vec, optimisation_method, end_criteria);
    
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            calibrated_g2[i*m+j][4]=swaption_helper_vec[i*m+j]->modelValue();
        }
    }
    std::cout << "Vasicek Params" << std::endl;
    std::cout <<vasicek->params()[0]<<std::endl;
    std::cout <<vasicek->params()[1]<<std::endl;
    std::cout <<vasicek->params()[2]<<std::endl;
    std::cout <<vasicek->params()[3]<<std::endl;
    std::cout << std::endl;
    
    writeCsv(calibrated_g2, "calibrated_g2.csv");
    
    
    Matrix calibrated_g2_term(300,5);
    
    t = 0; //already defined earlier
    for(int i=0;i<300;i++) {
        t+=0.1;
        calibrated_g2_term[i][0]=t;
        calibrated_g2_term[i][1]=the_handle->discount(t);
        calibrated_g2_term[i][2]=g2_ptr->discount(t);
        calibrated_g2_term[i][3]=hull_white->discount(t);
        calibrated_g2_term[i][4]=vasicek->discount(t);
    }
    
    writeCsv(calibrated_g2_term,"calibrated_g2_term.csv");
}

#endif /* CalibrateG2_hpp */
