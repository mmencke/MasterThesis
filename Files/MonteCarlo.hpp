//
//  McCva.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 10/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef McCva_hpp
#define McCva_hpp

#include "headers.hpp"


using namespace QuantLib;

Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap> > swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<G2> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery);

inline void mcCva() {
    
    
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
    
    std::vector<Date> the_dates;
    std::vector<Rate> the_rates;
    
    Matrix swapRatesData = readCsv("/Users/mmencke/Speciale_Data/swapRates.csv");
    
    
    Frequency fixedLegFrequency = Annual;
    Period fixedLegTenor = 1*Years;
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
    
    
    //std::cout <<yieldCurveHandle->forwardRate(1.81066,1.81066,Continuous, NoFrequency) << std::endl;
    
    /*******************************************************************/
    /* SHORT RATE MODELS                             */
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
    
    Matrix cdsSpreadBank = readCsv("/Users/mmencke/Speciale_Data/cdsSpreadBank.csv");
    
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
    
    
    
    
    Matrix cdsSpreadCpt = readCsv("/Users/mmencke/Speciale_Data/cdsSpreadCounterparty.csv");
    
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
    /* CREATING THE SWAPS                              */
    /*******************************************************************/
    
    
    Date startDate = calendar.advance(anchorDate,2, Days,ModifiedFollowing);
    Real nominal = 1000000.0;
    VanillaSwap::Type swapType = VanillaSwap::Payer;
    Rate tempFixedRate = 0.0;
    Spread swapSpread =0.0;
    std::vector<Period> exmpSwapTenors={5*Years, 10*Years, 15*Years, 20*Years, 25*Years, 30*Years};
    
    ext::shared_ptr<IborIndex> euriborIndex(new Euribor3M(yieldCurveHandle));
    
    euriborIndex->clearFixings(); //make sure that we don't have fixings from previous code
    

    std::vector<ext::shared_ptr<VanillaSwap>> atmSwapVec;
    std::vector<ext::shared_ptr<VanillaSwap>> otmSwapVec;
    std::vector<ext::shared_ptr<VanillaSwap>> itmSwapVec;
    
    Real rateAdj=0.0025; //25bp
    
    Matrix exmpSwaps(exmpSwapTenors.size(),7);
    
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
        
        VanillaSwap itmSwap=MakeVanillaSwap(exmpSwapTenors[i], euriborIndex,parRate-rateAdj)
        .withType(swapType)
        .withNominal(nominal)
        .withEffectiveDate(startDate)
        .withFloatingLegSpread(swapSpread)
        ;
        itmSwapVec.push_back(ext::make_shared<VanillaSwap>(itmSwap));
        
        VanillaSwap otmSwap=MakeVanillaSwap(exmpSwapTenors[i], euriborIndex,parRate+rateAdj)
        .withType(swapType)
        .withNominal(nominal)
        .withEffectiveDate(startDate)
        .withFloatingLegSpread(swapSpread)
        ;
        otmSwapVec.push_back(ext::make_shared<VanillaSwap>(otmSwap));
        
        exmpSwaps[i][0]=years(exmpSwapTenors[i]);
        exmpSwaps[i][1]=atmSwapVec[i]->fixedRate()*10000;
        exmpSwaps[i][2]=atmSwapVec[i]->NPV();
        exmpSwaps[i][3]=itmSwapVec[i]->fixedRate()*10000;
        exmpSwaps[i][4]=itmSwapVec[i]->NPV();
        exmpSwaps[i][5]=otmSwapVec[i]->fixedRate()*10000;
        exmpSwaps[i][6]=otmSwapVec[i]->NPV();

        //std::cout << exmpSwapTenors[i] << ": " <<atmSwapVec[i]->fixedRate() << ", " << itmSwapVec[i]->fixedRate() << ", " << otmSwapVec[i]->fixedRate() << std::endl;
        //std::cout << "   " <<atmSwapVec[i]->NPV() << ", " << itmSwapVec[i]->NPV() << ", " << otmSwapVec[i]->NPV() << std::endl;
    }
    
    writeCsv(exmpSwaps, "exmpSwaps.csv");
    
    
    Real cptRecovery = 0.4;
    Real bankRecovery = 0.4;
    
    unsigned long seed =1;
    
    
    Matrix rho(4,4,0);
    rho[0][0]=rho[1][1]=rho[2][2]=rho[3][3]=1;
    
    ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator(new MersenneTwisterUniformRng(seed));
    
    Size nSamples=10000;
    
    Real qNorm95=1.644853626951472;//call sprintf("%.15f",qnorm(0.95)) in R
    
    
    std::vector<ext::shared_ptr<VanillaSwap>> tempSvapVec;
    Matrix tempMatrix;

    /*******************************************************************/
    /* INDEPENDENCE                              */
    /*******************************************************************/
    

    start = std::chrono::steady_clock::now();
    Matrix cvaDvaInd(exmpSwapTenors.size(),13);

    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaInd[i][0]=years(exmpSwapTenors[i]);
        
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
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaInd[i][5]=tempMatrix[0][0];
        cvaDvaInd[i][6]= qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][7]=tempMatrix[0][1];
        cvaDvaInd[i][8]=qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaInd[i][9]=tempMatrix[0][0];
        cvaDvaInd[i][10]= qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][11]=tempMatrix[0][1];
        cvaDvaInd[i][12]=qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaInd,"cvaDvaInd.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;

    std::cout << "Independent: Run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    /*******************************************************************/
    /* WRONG WAY RISK                              */
    /*******************************************************************/

    rho[0][2]=rho[2][0]=-0.05;
    rho[1][2]=rho[2][1]=0.7505596;
    
    //printMatrix(rho);
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaWwr(exmpSwapTenors.size(),13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaWwr[i][0]=years(exmpSwapTenors[i]);
        
        Size nTimeSteps=52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][1]=tempMatrix[0][0];
        cvaDvaWwr[i][2]= qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][3]=tempMatrix[0][1];
        cvaDvaWwr[i][4]=qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][5]=tempMatrix[0][0];
        cvaDvaWwr[i][6]= qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][7]=tempMatrix[0][1];
        cvaDvaWwr[i][8]=qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][9]=tempMatrix[0][0];
        cvaDvaWwr[i][10]= qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][11]=tempMatrix[0][1];
        cvaDvaWwr[i][12]=qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaWwr,"cvaDvaWwr.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Wrong Way Risk: Run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    

    /*******************************************************************/
    /* RIGHT WAY RISK                              */
    /*******************************************************************/

    //rho[0][2]=rho[2][0]=-0.05;
    
    rho[1][2]=rho[2][1]=-0.648994;
    //printMatrix(rho);
    
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaRwr(exmpSwapTenors.size(),13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaRwr[i][0]=years(exmpSwapTenors[i]);
        
        Size nTimeSteps=52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][1]=tempMatrix[0][0];
        cvaDvaRwr[i][2]= qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][3]=tempMatrix[0][1];
        cvaDvaRwr[i][4]=qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][5]=tempMatrix[0][0];
        cvaDvaRwr[i][6]= qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][7]=tempMatrix[0][1];
        cvaDvaRwr[i][8]=qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix=mcCvaSwapPortfolio(tempSvapVec,uniformGenerator,nTimeSteps,nSamples,
                                      rho,g2Model, cirCptModel, cirBankModel,
                                      cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][9]=tempMatrix[0][0];
        cvaDvaRwr[i][10]= qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][11]=tempMatrix[0][1];
        cvaDvaRwr[i][12]=qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaRwr,"cvaDvaRwr.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Right Way Risk: Run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    //ext::shared_ptr<TreeSwaptionEngine> g2EnginePtr(new TreeSwaptionEngine(g2Model, 25));
    
/*
    ext::shared_ptr<G2SwaptionEngine> g2EnginePtr(new G2SwaptionEngine(g2Model, 4, 100));//±4 standard deviations wide with 100 intervals
    Handle<PricingEngine> g2EngineHandle(g2EnginePtr);
    
    ext::shared_ptr<PricingEngine> cvaAdjSwapEngine =
    ext::make_shared<CounterpartyAdjSwapEngine>(yieldCurveHandle,
                                                g2EngineHandle,
                                                defaultCurveHandleBank,
                                                bankRecovery
                                                );
    
    ext::shared_ptr<PricingEngine> bvaAdjSwapEngine =
    ext::make_shared<CounterpartyAdjSwapEngine>(
                                                yieldCurveHandle,
                                                g2EngineHandle,
                                                defaultCurveHandleCpt,
                                                cptRecovery,
                                                defaultCurveHandleBank,
                                                bankRecovery
                                                );
    
    ext::shared_ptr<VanillaSwap> riskySwapPtr(itmSwapVec[5]); //clone the swap to ptr
    
    riskySwapPtr->setPricingEngine(cvaAdjSwapEngine);
    
    Real noCrNpv=itmSwapVec[5]->NPV();
    
    Real cvaAdjNpv=riskySwapPtr->NPV();
    Real cva=noCrNpv-cvaAdjNpv;
    
    std::cout << "RF = " << noCrNpv << std::endl;
    
    std::cout << "CVA = " << cva << std::endl;
    

    riskySwapPtr->setPricingEngine(bvaAdjSwapEngine);
    
    Real bvaAdjNpv=riskySwapPtr->NPV();
    
    Real bva=bvaAdjNpv-noCrNpv;
    
    Real dva=bva+cva;
    
    std::cout << "DVA = " << dva << std::endl;
    std::cout << "BVA = " << bva << std::endl;
*/
    
    
    
    
    
    
    /*
    
    
    std::cout <<exmpSwap.fairRate() << std::endl;
    
    ext::shared_ptr<VanillaSwap> exmpSwapPtr=ext::make_shared<VanillaSwap>(exmpSwap);
    
    std::cout << exmpSwap.NPV() << std::endl;
    std::cout << exmpSwapPtr->NPV() << std::endl;
    
    std::vector<ext::shared_ptr<VanillaSwap> > swapPortfolio;
    
    swapPortfolio.push_back(exmpSwapPtr);
    
    std::cout <<swapPortfolio[0]->NPV()<<std::endl;

    
    Size nTimeSteps=10;//3 month steps: 4*10
    Size nSamples=100000;
    
    ext::shared_ptr<G2> interestRate(new G2(yieldCurveHandle));
    
    Array g2Params={0.0183149,0.0106535,0.0183039,0.0110041,-0.729877};
    
    interestRate->setParams(g2Params);
    
    
    Matrix cdsSpreads = readCsv("/Users/mmencke/Speciale_Data/table12-4.csv");
    
    
    //printMatrix(cds_spreads);
    
    // market
    Natural settlementDays = 3;
    Real marketRecovery = 0.4;
    std::vector<Period> tenors;
    for(int i=0;i<10;i++) {
        tenors.push_back((i+1) * Years);
    }
    
    
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > instruments;
    for (Size i = 0; i < tenors.size(); i++) {
        instruments.push_back(ext::shared_ptr<DefaultProbabilityHelper>(new SpreadCdsHelper(Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(cdsSpreads[i][0]))),
                                                                                            tenors[i], settlementDays, zeroCalendar, Quarterly, Following,
                                                                                            DateGeneration::CDS2015,Actual360(),
                                                                                            marketRecovery, yieldCurveHandle)));
    }
    // Bootstrap hazard rates
    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
    hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
                                                                            anchorDate, instruments, zeroDayCounter));
    
    RelinkableHandle<DefaultProbabilityTermStructure> defaultCurveHandle;
    defaultCurveHandle.linkTo(hazardRateStructure);

    std::vector<Date> hazardDates = hazardRateStructure->dates();
    //std::vector<Real> hazardRates = hazardRateStructure->data();
    
    std::vector<Real> discFactors;
    for(int i=0;i<hazardDates.size();i++) {
        discFactors.push_back(hazardRateStructure->survivalProbability(hazardDates[i]));
    }
    
    ext::shared_ptr<InterpolatedDiscountCurve<Linear>> equivalentYieldPtr(new InterpolatedDiscountCurve<Linear>(hazardDates,discFactors,zeroDayCounter));
    equivalentYieldPtr->enableExtrapolation();
    
    Handle<YieldTermStructure> equivalentYieldHandle(equivalentYieldPtr);
    
    ext::shared_ptr<ExtendedCoxIngersollRoss> ctptyIntensity(new ExtendedCoxIngersollRoss(equivalentYieldHandle));

    ext::shared_ptr<ExtendedCoxIngersollRoss> invstIntensity(new ExtendedCoxIngersollRoss(equivalentYieldHandle));
    
    Array extCirParams={0.0331913,0.780113,0.239354,0.0166051};

    ctptyIntensity->setParams(extCirParams);
    invstIntensity->setParams(extCirParams);
    
    Real ctptyRecovery = 0.4;
    Real invstRecovery = 0.4;

    unsigned long seed =1;

    
    Matrix rho(4,4,0);
    rho[0][0]=rho[1][1]=rho[2][2]=rho[3][3]=1;
 
    ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator(new MersenneTwisterUniformRng(seed));
    
    Matrix cvaDva=mcCvaSwapPortfolio(swapPortfolio,uniformGenerator,nTimeSteps,nSamples,
                                      rho,interestRate, ctptyIntensity, invstIntensity,
                                      ctptyRecovery, invstRecovery);
    
    
    std::cout << std::endl;
    printMatrix(cvaDva);
    std::cout << std::endl;

    
    
    //ext::shared_ptr<G2SwaptionEngine> swaptionEnginePtr(new G2SwaptionEngine(interestRate, 4, 100));//±4 standard deviations wide with 100 intervals
    ext::shared_ptr<TreeSwaptionEngine> swaptionEnginePtr(new TreeSwaptionEngine(interestRate, 25));
    Handle<PricingEngine> swaptionEngineHandle(swaptionEnginePtr);
    
    ext::shared_ptr<PricingEngine> cvaSwapEngine =
    ext::make_shared<CounterpartyAdjSwapEngine>(yieldCurveHandle,
                                                swaptionEngineHandle,
                                                defaultCurveHandle,
                                                ctptyRecovery
                                                );
    
    ext::shared_ptr<VanillaSwap> riskySwapPtr(new VanillaSwap(exmpSwap)); //clone the swap to ptr
    
    riskySwapPtr->setPricingEngine(cvaSwapEngine);
    
    Real rfNpv=exmpSwapPtr->NPV();
     
    Real cvaAdjusted=riskySwapPtr->NPV();
    
    std::cout << "CVA = " << -(cvaAdjusted-rfNpv) << std::endl;
    
    ext::shared_ptr<PricingEngine> cvaDvaSwapEngine =
    ext::make_shared<CounterpartyAdjSwapEngine>(
                                                yieldCurveHandle,
                                                swaptionEngineHandle,
                                                defaultCurveHandle,
                                                ctptyRecovery,
                                                defaultCurveHandle,
                                                invstRecovery
                                                );
    riskySwapPtr->setPricingEngine(cvaDvaSwapEngine);
    
    
    Real cvaDvaAdjusted=riskySwapPtr->NPV();
    
    std::cout << "DVA = " << cvaDvaAdjusted -cvaAdjusted << std::endl;
    
    
    
    std::cout << "Adjusted = " << cvaDvaAdjusted << ", " <<rfNpv-cvaDva[0][0]+cvaDva[0][1] << std::endl;*/
}


#endif /* McCva_hpp */
