//
//  Calibration.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 12/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef Calibration_hpp
#define Calibration_hpp

#include "headers.hpp"

inline void calibration() {
    
    Date anchorDate(04,May,2021);
    Settings::instance().evaluationDate() = anchorDate;
    
    /*******************************************************************/
    /* BOOTSTRAPPING YIELD CURVE                                       */
    /*******************************************************************/
    
    
    Actual360 floatingLegDayCounter;
    TARGET calendar;
    
    Matrix swapRatesData = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/swapRates.csv");
    
    
    Frequency fixedLegFrequency = Annual;
    Period fixedLegTenor = 1*Years;
    BusinessDayConvention fixedLegConvention = ModifiedFollowing;
    DayCounter fixedLegDayCounter = Thirty360(Thirty360::European);
    
    
    RelinkableHandle<YieldTermStructure> yieldCurveHandle;
    ext::shared_ptr<IborIndex> floatingLegIndex(new Euribor3M(yieldCurveHandle));
    
    std::vector<Period> swapRateTenors={1*Years,    2*Years,    3*Years,    4*Years,    5*Years,
        6*Years,    7*Years,    8*Years,    9*Years,    10*Years,    11*Years,    12*Years,
        15*Years,    20*Years,    25*Years,    30*Years,    40*Years,    50*Years};

    std::vector<ext::shared_ptr<RateHelper>> rateHelperVec;
    for(int i=0;i<swapRateTenors.size();i++) {
        ext::shared_ptr<RateHelper> temp(new SwapRateHelper(
                                                            Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(swapRatesData[0][i]/100))),//rates in percent
                                                            swapRateTenors[i],
                                                            calendar, fixedLegFrequency,
                                                            fixedLegConvention, fixedLegDayCounter,
                                                            floatingLegIndex));
        rateHelperVec.push_back(temp);
        
    }
    
    ext::shared_ptr<YieldTermStructure> yieldCurvePtr(new PiecewiseYieldCurve<Discount,Cubic>(anchorDate, rateHelperVec, floatingLegDayCounter));
    yieldCurvePtr->enableExtrapolation();
    yieldCurveHandle.linkTo(yieldCurvePtr);
    
    /*******************************************************************/
    /* SWAPTION CALIBRATION                                            */
    /*******************************************************************/
    
    Matrix swaptionVols = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/swaptionVols.csv");
    
    std::vector<ext::shared_ptr<SwaptionHelper>> swaptionHelperVec;
    
    std::vector<Period> swaptionExpiries={
        1*Months, 2*Months, 3*Months, 6*Months, 9*Months,
        1*Years, 18*Months, 2*Years, 3*Years, 4*Years,
        5*Years, 6*Years, 7*Years, 8*Years, 9*Years,
        10*Years, 15*Years, 20*Years, 25*Years, 30*Years
    };
    std::vector<Period> swaptionTenors={
        1*Years, 2*Years, 3*Years, 4*Years, 5*Years,
        6*Years, 7*Years, 8*Years, 9*Years, 10*Years,
        15*Years, 20*Years, 25*Years, 30*Years
    };

    for(int i=0;i<swaptionExpiries.size();i++) {
        for(int j=0;j<swaptionTenors.size();j++) {
            ext::shared_ptr<SwaptionHelper> temp(new SwaptionHelper(swaptionExpiries[i], swaptionTenors[j],
                                                                    Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(swaptionVols[i][j]/10000))),//vol in basispoints
                                                                    floatingLegIndex, fixedLegTenor, fixedLegDayCounter, floatingLegDayCounter,  yieldCurveHandle,
                                                                    BlackCalibrationHelper::ImpliedVolError, Null<Real>(), 1.0, Normal));

            swaptionHelperVec.push_back(temp);
        }
    }

    LevenbergMarquardt  optimisationMethod(1.0e-8,1.0e-8,1.0e-8);//epsfcn,xtol,gtol
    EndCriteria endCriteria(1000, 100, 1e-6, 1e-8, 1e-8);//maxIterations,  maxStationaryStateIterations,  rootEpsilon,  functionEpsilon,  gradientNormEpsilon
    
    ext::shared_ptr<G2> g2Model(new G2(yieldCurveHandle));
    ext::shared_ptr<PricingEngine> g2Engine(new G2SwaptionEngine(g2Model,4, 100));//±4 standard deviations wide with 100 intervals
    //ext::shared_ptr<PricingEngine> g2SwaptionEngine(new TreeSwaptionEngine(g2Model,25));
    //ext::shared_ptr<PricingEngine> g2SwaptionEngine(new FdG2SwaptionEngine(g2Model));
    
    std::vector<ext::shared_ptr<CalibrationHelper>> tempHelperVec;
    for(int i=0; i<swaptionHelperVec.size();i++) {
        if(i==7||i==8) {
            //do nothing
            //problems with 1Mo-8Yr and 1Mo-9Yr for G2,i==7||i==8
        } else {
            swaptionHelperVec[i]->setPricingEngine(g2Engine);
            tempHelperVec.push_back(swaptionHelperVec[i]);
        }
    }
    auto start = std::chrono::steady_clock::now();
    g2Model->calibrate(tempHelperVec, optimisationMethod, endCriteria);
    auto end = std::chrono::steady_clock::now();
    
    auto diff = end - start;
    std::cout << boost::format("%-35s %-8s %-8s") % "G2: Calibration run time =" % std::chrono::duration <double, std::ratio<60>> (diff).count() % "minutes" << std::endl;
    
    
    std::cout << "G2: List of Parameters" << std::endl;
    std::cout << "   " << g2Model->params()[0]<<std::endl;
    std::cout << "   " << g2Model->params()[1]<<std::endl;
    std::cout << "   " << g2Model->params()[2]<<std::endl;
    std::cout << "   " << g2Model->params()[3]<<std::endl;
    std::cout << "   " << g2Model->params()[4]<<std::endl;
    std::cout << std::endl;
    
    Array g2Params={0.00549236,0.00970193,0.00550213,0.00955238,-0.759051};
    g2Model->setParams(g2Params);
    
    ext::shared_ptr<HullWhite> hwModel(new HullWhite(yieldCurveHandle));
    ext::shared_ptr<PricingEngine> hwEngine(new JamshidianSwaptionEngine(hwModel,yieldCurveHandle));
    
    
    tempHelperVec.clear();
    for(int i=0; i<swaptionHelperVec.size();i++) {
        swaptionHelperVec[i]->setPricingEngine(hwEngine);
        tempHelperVec.push_back(swaptionHelperVec[i]);
    }
    start = std::chrono::steady_clock::now();
    hwModel->calibrate(tempHelperVec, optimisationMethod, endCriteria);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cout << boost::format("%-35s %-8s %-8s") % "Hull-White: Calibration run time =" % std::chrono::duration <double, std::ratio<60>> (diff).count() % "minutes" << std::endl;
    
    std::cout << "Hull-White: List of Parameters" << std::endl;
    std::cout << "   " << hwModel->params()[0]<<std::endl;
    std::cout << "   " << hwModel->params()[1]<<std::endl;
    std::cout << std::endl;
    
    Array hwParams={0.00531006,0.00673971};
    hwModel->setParams(hwParams);
    
    Matrix swaptionFit(swaptionHelperVec.size(),5);
    
    int k=0;
    for(int i=0;i<swaptionExpiries.size();i++) {
        for(int j=0;j<swaptionTenors.size();j++) {
            swaptionFit[k][0]=years(swaptionExpiries[i]);
            swaptionFit[k][1]=years(swaptionTenors[j]);
            
            swaptionFit[k][2]=swaptionHelperVec[k]->volatility().currentLink()->value()*10000;
            
            swaptionHelperVec[k]->setPricingEngine(g2Engine);
            if(k==7||k==8) {
                swaptionFit[k][3]=0;
            } else {
                //swaptionFit[k][3]=swaptionHelperVec[k]->modelValue();
                swaptionFit[k][3]=swaptionHelperVec[k]->impliedVolatility(swaptionHelperVec[k]->modelValue(),
                                                                          1e-08,//accuracy
                                                                          100, //max evals
                                                                          1e-08, //min vol
                                                                          1)*10000; //max vol
            }
            
            swaptionHelperVec[k]->setPricingEngine(hwEngine);
            swaptionFit[k][4]=swaptionHelperVec[k]->impliedVolatility(swaptionHelperVec[k]->modelValue(),
                                                                      1e-08,//accuracy
                                                                      100, //max evals
                                                                      1e-08, //min vol
                                                                      1)*10000; //max vol
            
            k++;
        }
    }
    
    writeCsv(swaptionFit,"swaptionFit.csv");
    
    
    Matrix calibratedZeroCurve(500,4);
    
    double t = 0; //already defined earlier
    for(int i=0;i<500;i++) {
        t+=0.1;
        calibratedZeroCurve[i][0]=t;
        calibratedZeroCurve[i][1]=(-log(yieldCurveHandle->discount(t))/t)*100;
        calibratedZeroCurve[i][2]=(-log(g2Model->discount(t))/t)*100;
        calibratedZeroCurve[i][3]=(-log(hwModel->discount(t))/t)*100;
    }
    
    writeCsv(calibratedZeroCurve,"calibratedZeroCurve.csv");
    
    /*******************************************************************/
    /* BOOTSTRAPPING HAZARD RATES                              */
    /*******************************************************************/
    
    Matrix bankCdsSpreads = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/bankCdsSpreads.csv");
    
    std::vector<Period> cdsTenors = {
        6*Months, 1*Years, 2*Years, 3*Years, 4*Years,
        5*Years, 7*Years, 10*Years, 15*Years, 20*Years,
        30*Years
    };
    
    Natural cdsSettlementDays = 3;
    Real cdsRecoveryRate = 0.4;
    Frequency cdsFrequency = Quarterly;
    BusinessDayConvention cdsPaymentConv =Following;
    DateGeneration::Rule cdsDateGen =DateGeneration::CDS2015;
    Actual360 cdsDayCounter;
    
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > cdsInstruments;
    for (Size i = 0; i < cdsTenors.size(); i++) {
        cdsInstruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(
                                                                           new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(bankCdsSpreads[0][i]))),
                                                                                               cdsTenors[i], cdsSettlementDays, calendar,
                                                                                               cdsFrequency, cdsPaymentConv,
                                                                                               cdsDateGen,cdsDayCounter,
                                                                                               cdsRecoveryRate, yieldCurveHandle)));
    }
    
    
    ext::shared_ptr<PiecewiseDefaultCurve<SurvivalProbability, Cubic>> bankDefaultCurvePtr(
                                                                    new PiecewiseDefaultCurve<SurvivalProbability, Cubic>(anchorDate, cdsInstruments, cdsDayCounter));
    
    std::vector<Date> defaultCurveDates = bankDefaultCurvePtr->dates();
    
    std::vector<Real> equivalentDiscFactors;
    
    for(int i=0;i<defaultCurveDates.size();i++) {
        equivalentDiscFactors.push_back(bankDefaultCurvePtr->survivalProbability(defaultCurveDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> bankEquivYieldCurvePtr(new InterpolatedDiscountCurve<Cubic>(defaultCurveDates,equivalentDiscFactors,cdsDayCounter));
    
    Handle<YieldTermStructure> bankEquivYieldCurveHandle(bankEquivYieldCurvePtr);
    
    // Parameters from [Brigo, Morini and Pallavicini (2013), p. 125]
    Real kappa=0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    bool withFellerConstraint = true;
    
    ext::shared_ptr<CoxIngersollRoss> bankCirModel(new ExtendedCoxIngersollRoss(bankEquivYieldCurveHandle,
                                                                                theta,  kappa, sigma,lambda0,
                                                                                withFellerConstraint));
    

    
    
    Matrix cptCdsSpreads = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/cptCdsSpreads.csv");
    
    cdsInstruments.clear();
    for (Size i = 0; i < cdsTenors.size(); i++) {
        cdsInstruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(
                                                                           new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cptCdsSpreads[0][i]))),
                                                                                               cdsTenors[i], cdsSettlementDays, calendar,
                                                                                               cdsFrequency, cdsPaymentConv,
                                                                                               cdsDateGen,cdsDayCounter,
                                                                                               cdsRecoveryRate, yieldCurveHandle)));
    }
    
    ext::shared_ptr<PiecewiseDefaultCurve<SurvivalProbability, Cubic>> cptDefaultCurvePtr(
                                                                         new PiecewiseDefaultCurve<SurvivalProbability, Cubic>(anchorDate, cdsInstruments, cdsDayCounter));
    
    
    equivalentDiscFactors.clear();
    
    for(int i=0;i<defaultCurveDates.size();i++) {
        equivalentDiscFactors.push_back(cptDefaultCurvePtr->survivalProbability(defaultCurveDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> cptEquivYieldCurvePtr(new InterpolatedDiscountCurve<Cubic>(defaultCurveDates,equivalentDiscFactors,cdsDayCounter));
    
    Handle<YieldTermStructure> cptEquivYieldCurveHandle(cptEquivYieldCurvePtr);
    
    ext::shared_ptr<CoxIngersollRoss> cptCirModel(new ExtendedCoxIngersollRoss(cptEquivYieldCurveHandle,
                                                                                theta,  kappa, sigma,lambda0,
                                                                                withFellerConstraint));
    
    Matrix calibratedDefaultCurve(300,5);
    
    t=0; //defined earlier
    for(int i=0;i<300;i++) { //30 years
        t+=0.1;
        calibratedDefaultCurve[i][0]=t;
        calibratedDefaultCurve[i][1]=bankDefaultCurvePtr->survivalProbability(t);
        calibratedDefaultCurve[i][2]=bankCirModel->discount(t);
        calibratedDefaultCurve[i][3]=cptDefaultCurvePtr->survivalProbability(t);
        calibratedDefaultCurve[i][4]=cptCirModel->discount(t);
    }
    writeCsv(calibratedDefaultCurve,"calibratedDefaultCurve.csv");

}


#endif /* Calibration_hpp */
