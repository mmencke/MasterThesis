//
//  MonteCarloSample.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 30/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef MonteCarloSample_hpp
#define MonteCarloSample_hpp

#include "headers.hpp"
#include "MonteCarlo.hpp"

inline void monteCarloSample() {

    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    
    Date anchorDate(04,May,2021);
    Settings::instance().evaluationDate() = anchorDate;
    
    /*******************************************************************/
    /* BOOTSTRAPPING YIELD CURVE                              */
    /*******************************************************************/
    
    
    Actual360 floatingLegDayCounter;
    TARGET calendar;
    
    Matrix swapRatesData = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/swapRates.csv");
    
    
    Frequency fixedLegFrequency = Annual;
    BusinessDayConvention fixedLegConvention = ModifiedFollowing;
    DayCounter fixedLegDayCounter = Thirty360(Thirty360::European);
    
    
    RelinkableHandle<YieldTermStructure> yieldCurveHandle;
    ext::shared_ptr<IborIndex> floatingLegIndex(new Euribor3M(yieldCurveHandle));
    
    std::vector<Period> swapRateTenors={1*Years,    2*Years,    3*Years,    4*Years,    5*Years,
        6*Years, 7*Years, 8*Years,    9*Years,    10*Years,    11*Years,    12*Years,
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
    
    ext::shared_ptr<YieldTermStructure> yieldCurvePtr(new PiecewiseYieldCurve<Discount,Linear>(anchorDate, rateHelperVec, floatingLegDayCounter));
    yieldCurvePtr->enableExtrapolation();
    yieldCurveHandle.linkTo(yieldCurvePtr);
    
    /*******************************************************************/
    /* SHORT RATE MODELS                                               */
    /*******************************************************************/
    
    ext::shared_ptr<G2> g2Model(new G2(yieldCurveHandle));
    ext::shared_ptr<PricingEngine> g2Engine(new G2SwaptionEngine(g2Model,4, 100));//±4 standard deviations wide with 100 intervals
    
    
    Array g2Params={0.00549236,0.00970193,0.00550213,0.00955238,-0.759051};
    g2Model->setParams(g2Params);
    
    ext::shared_ptr<HullWhite> hwModel(new HullWhite(yieldCurveHandle));
    ext::shared_ptr<PricingEngine> hwEngine(new JamshidianSwaptionEngine(hwModel,yieldCurveHandle));
    
    Array hwParams={0.00531006,0.00673971};
    hwModel->setParams(hwParams);
    
    /*******************************************************************/
    /* BOOTSTRAPPING HAZARD RATES                              */
    /*******************************************************************/
    
    Matrix cdsSpreadBank = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/bankCdsSpreads.csv");
    
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
                                                                           new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cdsSpreadBank[0][i]))),
                                                                                               cdsTenors[i], cdsSettlementDays, calendar,
                                                                                               cdsFrequency, cdsPaymentConv,
                                                                                               cdsDateGen,cdsDayCounter,
                                                                                               cdsRecoveryRate, yieldCurveHandle)));
    }
    
    
    ext::shared_ptr<PiecewiseDefaultCurve<SurvivalProbability, Cubic>> defaultCurvePtrBank(
                                                                                           new PiecewiseDefaultCurve<SurvivalProbability, Cubic>(anchorDate, cdsInstruments, cdsDayCounter));
    
    RelinkableHandle<DefaultProbabilityTermStructure> defaultCurveHandleBank;
    defaultCurveHandleBank.linkTo(defaultCurvePtrBank);
    
    std::vector<Date> defaultCurveDates = defaultCurvePtrBank->dates();
    
    std::vector<Real> equivalentDiscFactors;
    
    for(int i=0;i<defaultCurveDates.size();i++) {
        equivalentDiscFactors.push_back(defaultCurvePtrBank->survivalProbability(defaultCurveDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> equivalentPtrBank(new InterpolatedDiscountCurve<Cubic>(defaultCurveDates,equivalentDiscFactors,cdsDayCounter));
    
    Handle<YieldTermStructure> equivalentHandleBank(equivalentPtrBank);
    
    //from page 125 in Brigo, Morini and Pallavicini
    Real kappa=0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    bool withFellerConstraint = true;
    
    ext::shared_ptr<CoxIngersollRoss> cirBankModel(new ExtendedCoxIngersollRoss(equivalentHandleBank,
                                                                                theta,  kappa, sigma,lambda0,
                                                                                withFellerConstraint));
    
    
    
    
    Matrix cdsSpreadCpt = readCsv("/Users/mmencke/Documents/GitHub/MasterThesis/Data/CppInput/cptCdsSpreads.csv");
    
    cdsInstruments.clear();
    for (Size i = 0; i < cdsTenors.size(); i++) {
        cdsInstruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(
                                                                           new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cdsSpreadCpt[0][i]))),
                                                                                               cdsTenors[i], cdsSettlementDays, calendar,
                                                                                               cdsFrequency, cdsPaymentConv,
                                                                                               cdsDateGen,cdsDayCounter,
                                                                                               cdsRecoveryRate, yieldCurveHandle)));
    }
    
    ext::shared_ptr<PiecewiseDefaultCurve<SurvivalProbability, Cubic>> defaultCurvePtrCpt(
                                                                                          new PiecewiseDefaultCurve<SurvivalProbability, Cubic>(anchorDate, cdsInstruments, cdsDayCounter));
    
    RelinkableHandle<DefaultProbabilityTermStructure> defaultCurveHandleCpt;
    defaultCurveHandleCpt.linkTo(defaultCurvePtrCpt);
    
    equivalentDiscFactors.clear();
    
    for(int i=0;i<defaultCurveDates.size();i++) {
        equivalentDiscFactors.push_back(defaultCurvePtrCpt->survivalProbability(defaultCurveDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> equivalentPtrCpt(new InterpolatedDiscountCurve<Cubic>(defaultCurveDates,equivalentDiscFactors,cdsDayCounter));
    
    Handle<YieldTermStructure> equivalentHandleCpt(equivalentPtrCpt);
    
    ext::shared_ptr<CoxIngersollRoss> cirCptModel(new ExtendedCoxIngersollRoss(equivalentHandleCpt,
                                                                               theta,  kappa, sigma, lambda0,
                                                                               withFellerConstraint));
    
    /*******************************************************************/
    /* CREATING THE SWAPS                                              */
    /*******************************************************************/
    
    Date startDate = calendar.advance(anchorDate,2, Days,ModifiedFollowing);
    Real nominal = 1000000.0;
    VanillaSwap::Type swapType = VanillaSwap::Payer;
    Rate tempFixedRate = 0.0;
    Spread swapSpread =0.0;
    
    
    std::vector<Period> exmpSwapTenors={30*Years};
    
    
    ext::shared_ptr<IborIndex> euriborIndex(new Euribor3M(yieldCurveHandle));
    
    euriborIndex->clearFixings(); //make sure that we don't have fixings from previous code
    
    
    std::vector<ext::shared_ptr<VanillaSwap>> atmSwapVec;
    std::vector<ext::shared_ptr<VanillaSwap>> otmSwapVec;
    std::vector<ext::shared_ptr<VanillaSwap>> itmSwapVec;
    
    Real rateAdj=0.0025; //25bp
    
    Matrix exmpSwaps(exmpSwapTenors.size(),3);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        VanillaSwap tempSwap=MakeVanillaSwap(exmpSwapTenors[i], euriborIndex,tempFixedRate)
        .withType(swapType)
        .withNominal(nominal)
        .withEffectiveDate(startDate)
        .withFloatingLegSpread(swapSpread)
        ;
        Real parRate = tempSwap.fairRate();
        
        VanillaSwap atmSwap=MakeVanillaSwap(exmpSwapTenors[i], euriborIndex,parRate)
        .withType(swapType)
        .withNominal(nominal)
        .withEffectiveDate(startDate)
        .withFloatingLegSpread(swapSpread)
        ;
        atmSwapVec.push_back(ext::make_shared<VanillaSwap>(atmSwap));
        
        exmpSwaps[i][0]=years(exmpSwapTenors[i]);
        exmpSwaps[i][1]=atmSwapVec[i]->fixedRate()*10000;
        exmpSwaps[i][2]=atmSwapVec[i]->NPV();
    }
    
    std::cout << std::endl;
    printMatrix(exmpSwaps);
    std::cout << std::endl;
    
    Real cptRecovery = 0.4;
    Real bankRecovery = 0.4;
    
    unsigned long seed =1;
    
    Matrix rho(4,4,0);
    rho[0][0]=rho[1][1]=rho[2][2]=rho[3][3]=1;
    
    ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator(new MersenneTwisterUniformRng(seed));
    
    Size nSamples=10000;
    
    Real qNorm95=1.959963984540054;//Source: sprintf("%.15f",qnorm(0.975)) in R
    
    
    std::vector<ext::shared_ptr<VanillaSwap>> tempSvapVec;
    Matrix tempMatrix;
    
    /*******************************************************************/
    /* INDEPENDENCE                                                    */
    /*******************************************************************/
    
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaInd(exmpSwapTenors.size(),5);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaInd[i][0]=years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps=52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaInd[i][1]=tempMatrix[0][0];
        cvaDvaInd[i][2]= qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][3]=tempMatrix[0][1];
        cvaDvaInd[i][4]=qNorm95*tempMatrix[1][1];
        
    }
    
    std::cout << std::endl;
    printMatrix(cvaDvaInd);
    std::cout << std::endl;
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Monte Carlo Sample: Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
}

#endif /* MonteCarloSample_hpp */
