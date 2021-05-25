//
//  CalibrateCir.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 28/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef CalibrateCir_hpp
#define CalibrateCir_hpp

#include "headers.hpp"

inline void calibrateCir() {
    Date anchorDate(26,May,2009);
    Settings::instance().evaluationDate() = anchorDate;
    
    Actual360 cdsDayCount;
    Actual360 zeroDayCount;
    TARGET calendar;
    
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    
    
    /*******************************************************************/
    /* CREATE YIELD CURVE OBJECT                                       */
    /*******************************************************************/
    
    std::vector<Date> zeroDates;
    std::vector<Rate> zeroRates;
    
    
    // Table 12.1 in [Brigo, Morini and Pallavicini (2013), p. 286]
    Matrix zeroTermStructure = readCsv("/Users/mmencke/Speciale_Data/table12-1.csv");
    
    for(int i=0;i<zeroTermStructure.rows();i++) {
        zeroDates.push_back(Date(zeroTermStructure[i][0]));
        zeroRates.push_back(zeroTermStructure[i][1]);
    }
    
    ext::shared_ptr<YieldTermStructure> yieldCurvePtr(new InterpolatedZeroCurve<Cubic>(zeroDates, zeroRates, zeroDayCount, calendar));
    
    RelinkableHandle<YieldTermStructure> yieldCurveHandle;
    yieldCurveHandle.linkTo(yieldCurvePtr);
    
    
    /*******************************************************************/
    /* BOOTSTRAP DEFAULT PROBABILITY CURVE                             */
    /*******************************************************************/
    
    // Table 12.4 in [Brigo, Morini and Pallavicini (2013), p. 289]
    Matrix cdsSpreads = readCsv("/Users/mmencke/Speciale_Data/table12-4.csv");
    
    Natural settlementDays = 3;
    Real recoveryRate = 0.4;
    Frequency cdsFrequency = Quarterly;
    DateGeneration::Rule cdsDateGen=DateGeneration::CDS2015;
    BusinessDayConvention cdsConvention=Following;
    
    std::vector<Period> tenors;
    for(int i=0;i<10;i++) {
        tenors.push_back((i+1) * Years);
    }
    
    
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > instruments;
    for (Size i = 0; i < tenors.size(); i++) {
        instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cdsSpreads[i][0]))),
                                                                                            tenors[i], settlementDays, calendar, cdsFrequency, cdsConvention,
                                                                                            cdsDateGen,cdsDayCount,
                                                                                            recoveryRate, yieldCurveHandle)));
    }

    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
    hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
                                                                            anchorDate, instruments, cdsDayCount));
    
    RelinkableHandle<DefaultProbabilityTermStructure> defaultCurveHandle;
    defaultCurveHandle.linkTo(hazardRateStructure);
    
    
    /*******************************************************************/
    /* CREATE CDSO CALIBRATION HELPER                                  */
    /*******************************************************************/
    
    Matrix cdsoVol = readCsv("/Users/mmencke/Speciale_Data/table12-7.csv");
    
    std::vector<ext::shared_ptr<CdsOptionHelper>> cdsoHelperVec;
    
    std::cout << "CDS Option Volatility: " << std::endl << std::endl;
    
    int k=0;
    std::cout << boost::format("%-6s %-3s %-6s %-3s %6s") % "Expiry" % " | " % "Length" % " | " % "Vol." << std::endl;
    std::cout << "------" << "---" << "------" << "---" << "------" << std::endl;
    
    for(int i=1; i<=10;i++) {
        for(int j=10;j>=1+i;j--) {
            std::cout << boost::format("%-6s %-3s %-6s %-3s %6s") % i % " | " % (10-(j-1)) % " | " % cdsoVol[k][0] << std::endl;
            
            ext::shared_ptr<CdsOptionHelper> tempHelperPtr(new CdsOptionHelper(i*Years,(10-(j-1))*Years,
                                                                               Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cdsoVol[k][0]))), recoveryRate, defaultCurveHandle, yieldCurveHandle));
            
            
            cdsoHelperVec.push_back(tempHelperPtr);
            
            k++; //counter for the matrix
        }
    }
    std::cout << std::endl;
    
    /*******************************************************************/
    /* CALIBRATE CIR                                                   */
    /*******************************************************************/
    
    
    //Parameters from [Brigo, Morini and Pallavicini (2013), p. 125]
    //Only used initially. Parameters will be changed in calibration
    Real kappa=0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    BigNatural seed = 1;
    Size nSamples = 10000;
    
    ext::shared_ptr<CoxIngersollRoss> cirModel(new CoxIngersollRoss(lambda0, theta, kappa, sigma));
    ext::shared_ptr<PricingEngine> cirCdsoEngine(new McCirCdsOptionEngine(cirModel, nSamples, seed, recoveryRate, yieldCurveHandle));
    
    LevenbergMarquardt  optMethod(1.0e-8,1.0e-8,1.0e-8);//epsfcn,xtol,gtol
    EndCriteria endCriteria(1000, 100, 1e-6, 1e-8, 1e-8);//maxIterations,  maxStationaryStateIterations,  rootEpsilon,  functionEpsilon,  gradientNormEpsilon
    
    std::vector<ext::shared_ptr<CalibrationHelper>> tempHelperVec;
    for(int k=0; k<cdsoHelperVec.size();k++){
        cdsoHelperVec[k]->setPricingEngine(cirCdsoEngine);
        tempHelperVec.push_back(cdsoHelperVec[k]);
    }
    
    start = std::chrono::steady_clock::now();
    cirModel->calibrate(tempHelperVec, optMethod, endCriteria);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    
    std::cout << "CIR: Run time = " << std::chrono::duration <double, std::ratio<60>>(diff).count() << " minutes" << std::endl;
    
    std::cout << "CIR: Parameters = " << std::endl;
    std::cout << "   " << cirModel->params()[0] << std::endl;
    std::cout << "   " << cirModel->params()[1] << std::endl;
    std::cout << "   " << cirModel->params()[2] << std::endl;
    std::cout << "   " << cirModel->params()[3] << std::endl;
    std::cout << std::endl;
    
    Array cirParams={0.0383418, 0.746503, 0.144219, 0.0397639};
    cirModel->setParams(cirParams);

    Matrix calibratedCir(cdsoHelperVec.size(),5);
    
    
    k=0;
    
    for(int i=1; i<=10;i++) {
        for(int j=10;j>=1+i;j--) {
            calibratedCir[k][0]=i;
            calibratedCir[k][1]=10-(j-1);
            
            k++; //counter for the matrix
        }
    }
    
    
    for(int k=0; k<cdsoHelperVec.size();k++){
        calibratedCir[k][2]=cdsoHelperVec[k]->marketValue();
        calibratedCir[k][3]=cdsoHelperVec[k]->modelValue();
    }
    
    
    /*******************************************************************/
    /* TRANSFORMING THE DEFAULT CURVE INTO AN EQUIVALENT YIELD CURVE   */
    /*******************************************************************/
    
    std::vector<Date> defaultDates = hazardRateStructure->dates();
    std::vector<Real> defaultRates = hazardRateStructure->data();
    
    std::vector<Real> equivDiscFactors;
    for(int i=0;i<defaultDates.size();i++) {
        equivDiscFactors.push_back(hazardRateStructure->survivalProbability(defaultDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Linear>> equivYieldCurvePtr(new InterpolatedDiscountCurve<Linear>(defaultDates,equivDiscFactors,cdsDayCount));
    
    Handle<YieldTermStructure> equivYieldCurveHandle(equivYieldCurvePtr);
    
    std::cout << "Check that the discount factor from the equivalent yield curve is equal to the survival probability from the default curve:" << std::endl << std::endl;
    std::cout << boost::format("%-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s") % "Date 1" % " | " % "Date 2" % " | " % "Rate 1" % " | " % "Rate 2" % " | " % "ZCB 1" % " | " % "ZCB 2" << std::endl;
    std::cout << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << "----" << "-----------" << std::endl;
    for(int i=0;i<defaultDates.size();i++) {
        Real equivZeroRate =equivYieldCurvePtr->zeroRate(defaultDates[i],cdsDayCount,Continuous);
        std::cout << boost::format("%-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s %-3s %-10s") % defaultDates[i] % " | " % equivYieldCurvePtr->dates()[i] % " | " % defaultRates[i] % " | " % equivZeroRate % " | " % hazardRateStructure->survivalProbability(defaultDates[i]) % " | " % equivYieldCurvePtr->discount(defaultDates[i]) << std::endl;
    }
    std::cout << std::endl;
    
    
    /*******************************************************************/
    /* CALIBRATE CIR++                                                 */
    /*******************************************************************/
    
    // Here we start calibration in the calibrated parameters from CIR
    
    ext::shared_ptr<ExtendedCoxIngersollRoss> extCirModel(new ExtendedCoxIngersollRoss(equivYieldCurveHandle,cirParams[0],cirParams[1],cirParams[2],cirParams[3]));
    
    ext::shared_ptr<PricingEngine> extCirCdsoEngine(new McCirCdsOptionEngine(extCirModel, nSamples, seed, recoveryRate, yieldCurveHandle));
    
    
    tempHelperVec.clear();
     for(int k=0; k<cdsoHelperVec.size();k++){
        cdsoHelperVec[k]->setPricingEngine(extCirCdsoEngine);
        tempHelperVec.push_back(cdsoHelperVec[k]);
    }
    
    start = std::chrono::steady_clock::now();
    extCirModel->calibrate(tempHelperVec, optMethod, endCriteria);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    
    std::cout << "CIR++: Run time = " << std::chrono::duration <double, std::ratio<60>>(diff).count() << " minutes" << std::endl;
    
    
    std::cout << "CIR++: Parameters = " << std::endl;
    std::cout << "   " << extCirModel->params()[0] << std::endl;
    std::cout << "   " << extCirModel->params()[1] << std::endl;
    std::cout << "   " << extCirModel->params()[2] << std::endl;
    std::cout << "   " << extCirModel->params()[3] << std::endl;
    std::cout << std::endl;
    
    
    
    Array extCirParams={0.033287, 0.782708, 0.239351, 0.0166414};
    extCirModel->setParams(extCirParams);
    
    for(int k=0; k<cdsoHelperVec.size();k++){
        calibratedCir[k][4]=cdsoHelperVec[k]->modelValue();
    }
    
    writeCsv(calibratedCir,"calibratedCir.csv");
    
    
    /*******************************************************************/
    /* CHECK TERM STRUCTURE                                            */
    /*******************************************************************/
    
    Matrix calibratedCirTerm(100,4);
    
    Time t = 0;
    for(int i=0;i<100;i++) {
        t+=0.1;
        calibratedCirTerm[i][0]=t;
        calibratedCirTerm[i][1]=hazardRateStructure->survivalProbability(t);
        calibratedCirTerm[i][2]=cirModel->discount(t);
        calibratedCirTerm[i][3]=extCirModel->discount(t);
    }
    
    writeCsv(calibratedCirTerm,"calibratedCirTerm.csv");
    
    
}


#endif /* CalibrateCir_hpp */
