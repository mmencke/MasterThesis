//
//  calibrateCir.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 28/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef calibrateCir_hpp
#define calibrateCir_hpp

#include "headers.hpp"

inline void calibrateCir() {
    Date anchor_date(26,May,2009);
    Settings::instance().evaluationDate() = anchor_date;
    
    Actual360 the_day_count;
    Actual360 zero_day_count;
    TARGET the_calendar;
    
    std::vector<Date> the_dates;
    std::vector<Rate> the_rates;
    
    Matrix term_structure = readCsv("/Users/mmencke/Speciale_Data/table12-1.csv");
    
    for(int i=0;i<term_structure.rows();i++) {
        the_dates.push_back(Date(term_structure[i][0]));
        the_rates.push_back(term_structure[i][1]);
    }
    
    ext::shared_ptr<YieldTermStructure> yield_curve_ptr(new InterpolatedZeroCurve<Cubic>(the_dates,the_rates, zero_day_count,the_calendar));
    
    RelinkableHandle<YieldTermStructure> yield_curve_handle;
    yield_curve_handle.linkTo(yield_curve_ptr);
    
    Matrix cds_spreads = readCsv("/Users/mmencke/Speciale_Data/table12-4.csv");
    
    
    //printMatrix(cds_spreads);
    
    // market
    Natural settlementDays = 3;
    Real recovery_rate = 0.4;
    std::vector<Period> tenors;
    for(int i=0;i<10;i++) {
        tenors.push_back((i+1) * Years);
    }
    
    
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > instruments;
    for (Size i = 0; i < tenors.size(); i++) {
        instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cds_spreads[i][0]))),
                                                                                            tenors[i], settlementDays, the_calendar, Quarterly, Following,
                                                                                            DateGeneration::CDS2015,Actual360(),
                                                                                            recovery_rate, yield_curve_handle)));
    }
    // Bootstrap hazard rates
    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
    hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
                                                                            anchor_date, instruments, the_day_count));
    
    RelinkableHandle<DefaultProbabilityTermStructure> default_curve_handle;
    default_curve_handle.linkTo(hazardRateStructure);
    
    
    //from page 125 in Brigo, Morini and Pallavicini
    Real kappa=0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    BigNatural seed = 1;
    Size n_samples = 10000;
    
    ext::shared_ptr<CoxIngersollRoss> cir_model(new CoxIngersollRoss(lambda0, theta, kappa, sigma));
    ext::shared_ptr<PricingEngine> cdsOptionEngine(new McCirCdsOptionEngine(cir_model, n_samples, seed, recovery_rate, yield_curve_handle));
    
    //CdsOptionHelper test(1*Years,1*Years, vol_handle, 0.4, default_curve_handle, yield_curve_handle);
    
    
     Matrix csdo_vol = readCsv("/Users/mmencke/Speciale_Data/table12-7.csv");
    
    std::vector<ext::shared_ptr<CdsOptionHelper>> cds_option_helper_vec;
    
    int k=0;
    std::cout << std::endl;
    std::cout << boost::format("%-6s %-3s %-6s %-3s %6s") % "Expiry" % " | " % "Length" % " | " % "Vol." << std::endl;
    std::cout << "------" << "---" << "-------" << "---" << "-------" << std::endl;
    
    for(int i=1; i<=10;i++) {
        for(int j=10;j>=1+i;j--) {
            std::cout << boost::format("%-6s %-3s %-6s %-3s %6s") % i % " | " % (10-(j-1)) % " | " % csdo_vol[k][0] << std::endl;
            
            ext::shared_ptr<CdsOptionHelper> temp_helper_ptr(new CdsOptionHelper(i*Years,(10-(j-1))*Years,
                                                                               Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(csdo_vol[k][0]))), recovery_rate, default_curve_handle, yield_curve_handle));
            
            temp_helper_ptr->setPricingEngine(cdsOptionEngine);
            
            cds_option_helper_vec.push_back(temp_helper_ptr);
            
            k++; //counter for the matrix
        }
    }
    std::cout << std::endl;
    
    LevenbergMarquardt  optimisation_method(1.0e-8,1.0e-8,1.0e-8);//epsfcn,xtol,gtol
    EndCriteria end_criteria(1000, 100, 1e-6, 1e-8, 1e-8);//maxIterations,  maxStationaryStateIterations,  rootEpsilon,  functionEpsilon,  gradientNormEpsilon
    
    std::vector<ext::shared_ptr<CalibrationHelper>> temp;
    for(int k=0; k<cds_option_helper_vec.size();k++){
        temp.push_back(cds_option_helper_vec[k]);
    }
    
    //cir_model->calibrate(temp, optimisation_method, end_criteria);
    
    std::cout << "CIR Parameters" << std::endl;
    std::cout << cir_model->params()[0] << std::endl;
    std::cout << cir_model->params()[1] << std::endl;
    std::cout << cir_model->params()[2] << std::endl;
    std::cout << cir_model->params()[3] << std::endl;
    std::cout << std::endl;
    
    Array cir_params={0.0384301,0.745443,0.144208,0.0400496};
    cir_model->setParams(cir_params);
    
    int n=cds_option_helper_vec.size();
    
    Matrix calibrated_cir(n,5);
    
    
    k=0;
    
    for(int i=1; i<=10;i++) {
        for(int j=10;j>=1+i;j--) {
            calibrated_cir[k][0]=i;
            calibrated_cir[k][1]=10-(j-1);
            
            k++; //counter for the matrix
        }
    }
    
    
    for(int k=0; k<cds_option_helper_vec.size();k++){
        calibrated_cir[k][2]=cds_option_helper_vec[k]->marketValue();
        calibrated_cir[k][3]=cds_option_helper_vec[k]->modelValue();
    }
    
    //Table 12.5
    /*the_params[0]=0.03;
    the_params[1]=0.5;
    the_params[2]=0.5;
    the_params[3]=0.05;
    
    cir_model->setParams(the_params);
    

    for(int k=0; k<cds_option_helper_vec.size();k++) {
        calibrated_cir[k][4]=cds_option_helper_vec[k]->modelValue();
    }
    */
    
    
    std::vector<Date> curve_dates = hazardRateStructure->dates();
    std::vector<Real> curve_rates = hazardRateStructure->data();
    
    std::vector<Real> disc_factors;
    for(int i=0;i<curve_dates.size();i++) {
        disc_factors.push_back(hazardRateStructure->survivalProbability(curve_dates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Linear>> equivalent_yield_curve(new InterpolatedDiscountCurve<Linear>(curve_dates,disc_factors,the_day_count));
    
    Handle<YieldTermStructure> equivalent_handle(equivalent_yield_curve);
    
    std::cout << std::endl;
    std::cout << boost::format("%-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s") % "Date 1" % " | " % "Date 2" % " | " % "Rate 1" % " | " % "Rate 2" % " | " % "ZCB 1" % " | " % "ZCB 2" << std::endl;
    std::cout << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << std::endl;
    for(int i=0;i<curve_dates.size();i++) {
        Real zero_rate =equivalent_yield_curve->zeroRate(curve_dates[i],the_day_count,Continuous);
        std::cout << boost::format("%-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s") % curve_dates[i] % " | " % equivalent_yield_curve->dates()[i] % " | " % curve_rates[i] % " | " % zero_rate % " | " % hazardRateStructure->survivalProbability(curve_dates[i]) % " | " % equivalent_yield_curve->discount(curve_dates[i]) << std::endl;
    }
    std::cout << std::endl;
    
    ext::shared_ptr<ExtendedCoxIngersollRoss> extCir(new ExtendedCoxIngersollRoss(equivalent_handle,cir_params[0],cir_params[1],cir_params[2],cir_params[3]));
    
    ext::shared_ptr<PricingEngine> cdsOptionEngineExt(new McCirCdsOptionEngine(extCir, n_samples, seed, recovery_rate, yield_curve_handle));
    
    
    
     for(int k=0; k<cds_option_helper_vec.size();k++){
        cds_option_helper_vec[k]->setPricingEngine(cdsOptionEngineExt);
    }
    

    for(int k=0; k<cds_option_helper_vec.size();k++){
        temp.push_back(cds_option_helper_vec[k]);
    }
    
    
    //extCir->calibrate(temp, optimisation_method, end_criteria);
    
    Array extCirParams={0.0331913,0.780113,0.239354,0.0166051};
    
    std::cout << "CIR++ Parameters" << std::endl;
    std::cout << extCir->params()[0] << std::endl;
    std::cout << extCir->params()[1] << std::endl;
    std::cout << extCir->params()[2] << std::endl;
    std::cout << extCir->params()[3] << std::endl;
    std::cout << std::endl;
    
    extCir->setParams(extCirParams);
    //extCir->setParams(cir_params);
    
    for(int k=0; k<cds_option_helper_vec.size();k++){
        calibrated_cir[k][4]=cds_option_helper_vec[k]->modelValue();
    }
    
    writeCsv(calibrated_cir,"calibrated_cir.csv");
    
    Matrix calibrated_cir_term(100,4);
    
    Time t = 0;
    for(int i=0;i<100;i++) {
        t+=0.1;
        calibrated_cir_term[i][0]=t;
        calibrated_cir_term[i][1]=hazardRateStructure->survivalProbability(t);
        //calibrated_cir_term[i][1]=equivalent_yield_curve->discount(t);
        calibrated_cir_term[i][2]=cir_model->discount(t);
        calibrated_cir_term[i][3]=extCir->discount(t);
    }
    
    writeCsv(calibrated_cir_term,"calibrated_cir_term.csv");
    
    
}


#endif /* calibrateCir_hpp */
