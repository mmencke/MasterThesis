//
//  MonteCarlo.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 10/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef MonteCarlo_hpp
#define MonteCarlo_hpp

#include "headers.hpp"


using namespace QuantLib;

Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap> > swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<G2> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery);

Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap> > swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, ext::shared_ptr<SobolRsg>, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<G2> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery);


Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap> > swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<HullWhite> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
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
    
    
    //std::cout <<yieldCurveHandle->forwardRate(1.81066,1.81066,Continuous, NoFrequency) << std::endl;
    
    /*******************************************************************/
    /* SHORT RATE MODELS                                               */
    /*******************************************************************/
    
    ext::shared_ptr<G2> g2Model(new G2(yieldCurveHandle));
 

    Array g2Params={0.00549236,0.00970193,0.00550213,0.00955238,-0.759051};
    g2Model->setParams(g2Params);
    
    ext::shared_ptr<HullWhite> hwModel(new HullWhite(yieldCurveHandle));
    
    Array hwParams={0.00531006,0.00673971};
    hwModel->setParams(hwParams);
    
    
    
    /*******************************************************************/
    /* BOOTSTRAPPING HAZARD RATES                                      */
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
    
    // Parameters from [Brigo, Morini and Pallavicini (2013), p. 125]
    Real kappa = 0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    bool withFellerConstraint = true;
    
    ext::shared_ptr<CoxIngersollRoss> cirBankModel(new ExtendedCoxIngersollRoss(equivalentHandleBank,
                                                                                theta,  kappa, sigma, lambda0,
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
    }
    
    writeCsv(exmpSwaps, "exmpSwaps.csv");
    
    
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
    /* INDEPENDENCE: G2++ PSEUDO-RANDOM                                */
    /*******************************************************************/

    start = std::chrono::steady_clock::now();
    Matrix cvaDvaInd(exmpSwapTenors.size(), 13);

    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaInd[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaInd[i][1] = tempMatrix[0][0];
        cvaDvaInd[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][3] = tempMatrix[0][1];
        cvaDvaInd[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaInd[i][5] = tempMatrix[0][0];
        cvaDvaInd[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][7] = tempMatrix[0][1];
        cvaDvaInd[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaInd[i][9] = tempMatrix[0][0];
        cvaDvaInd[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaInd[i][11] = tempMatrix[0][1];
        cvaDvaInd[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaInd, "cvaDvaInd.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;

    std::cout << "Independent (G2++ Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
 
    /*******************************************************************/
    /* INDEPENDENCE: G2++ QUASI-RANDOM                                 */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaIndSobol(exmpSwapTenors.size(), 13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaIndSobol[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        ext::shared_ptr<SobolRsg> sobolRsg(new SobolRsg(4*nTimeSteps,seed));
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndSobol[i][1] = tempMatrix[0][0];
        cvaDvaIndSobol[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaIndSobol[i][3] = tempMatrix[0][1];
        cvaDvaIndSobol[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndSobol[i][5] = tempMatrix[0][0];
        cvaDvaIndSobol[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaIndSobol[i][7] = tempMatrix[0][1];
        cvaDvaIndSobol[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndSobol[i][9] = tempMatrix[0][0];
        cvaDvaIndSobol[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaIndSobol[i][11] = tempMatrix[0][1];
        cvaDvaIndSobol[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaIndSobol, "cvaDvaIndSobol.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Independent (G2++ Quasi-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
    
    /*******************************************************************/
    /* INDEPENDENCE: HULL-WHITE PSEUDO-RANDOM                          */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaIndHw(exmpSwapTenors.size(), 13);
    
    Matrix rhoHw(3,3,0);
    rhoHw[0][0]=rhoHw[1][1]=rhoHw[2][2]=1;
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaIndHw[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndHw[i][1] = tempMatrix[0][0];
        cvaDvaIndHw[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaIndHw[i][3] = tempMatrix[0][1];
        cvaDvaIndHw[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndHw[i][5] = tempMatrix[0][0];
        cvaDvaIndHw[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaIndHw[i][7] = tempMatrix[0][1];
        cvaDvaIndHw[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaIndHw[i][9] = tempMatrix[0][0];
        cvaDvaIndHw[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaIndHw[i][11] = tempMatrix[0][1];
        cvaDvaIndHw[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaIndHw, "cvaDvaIndHw.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Independent (Hull-White Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
    
    /*******************************************************************/
    /* WRONG WAY RISK: G2++ PSEUDO-RANDOM                              */
    /*******************************************************************/

    rho[0][2]=rho[2][0]=-0.05;
    rho[1][2]=rho[2][1]=0.7505596;

    start = std::chrono::steady_clock::now();
    Matrix cvaDvaWwr(exmpSwapTenors.size(), 13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaWwr[i][0] = years(exmpSwapTenors[i]);
        
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][1] = tempMatrix[0][0];
        cvaDvaWwr[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][3] = tempMatrix[0][1];
        cvaDvaWwr[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][5] = tempMatrix[0][0];
        cvaDvaWwr[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][7] = tempMatrix[0][1];
        cvaDvaWwr[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwr[i][9] = tempMatrix[0][0];
        cvaDvaWwr[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaWwr[i][11] = tempMatrix[0][1];
        cvaDvaWwr[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaWwr, "cvaDvaWwr.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Wrong Way Risk (G2++ Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
 
    /*******************************************************************/
    /* WRONG WAY RISK: G2++ QUASI-RANDOM                               */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaWwrSobol(exmpSwapTenors.size(), 13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaWwrSobol[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        ext::shared_ptr<SobolRsg> sobolRsg(new SobolRsg(4*nTimeSteps, seed));
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrSobol[i][1] = tempMatrix[0][0];
        cvaDvaWwrSobol[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrSobol[i][3] = tempMatrix[0][1];
        cvaDvaWwrSobol[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrSobol[i][5] = tempMatrix[0][0];
        cvaDvaWwrSobol[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrSobol[i][7] = tempMatrix[0][1];
        cvaDvaWwrSobol[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho,g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrSobol[i][9] = tempMatrix[0][0];
        cvaDvaWwrSobol[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrSobol[i][11] = tempMatrix[0][1];
        cvaDvaWwrSobol[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaWwrSobol, "cvaDvaWwrSobol.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Wrong Way Risk (G2++ Quasi-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
    
    /*******************************************************************/
    /* WRONG WAY RISK: HULL-WHITE PSEUDO-RANDOM                        */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaWwrHw(exmpSwapTenors.size(),13);

    rhoHw[0][1]=rhoHw[1][0]=0.99;
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaWwrHw[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrHw[i][1] = tempMatrix[0][0];
        cvaDvaWwrHw[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrHw[i][3] = tempMatrix[0][1];
        cvaDvaWwrHw[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrHw[i][5] = tempMatrix[0][0];
        cvaDvaWwrHw[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrHw[i][7] = tempMatrix[0][1];
        cvaDvaWwrHw[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaWwrHw[i][9] = tempMatrix[0][0];
        cvaDvaWwrHw[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaWwrHw[i][11] = tempMatrix[0][1];
        cvaDvaWwrHw[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaWwrHw, "cvaDvaWwrHw.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Wrong Way Risk (Hull-White Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;

    /*******************************************************************/
    /* RIGHT WAY RISK: G2++ PSEUDO-RANDOM                              */
    /*******************************************************************/

    rho[1][2]=rho[2][1]=-0.648994;
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaRwr(exmpSwapTenors.size(),13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaRwr[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][1] = tempMatrix[0][0];
        cvaDvaRwr[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][3] = tempMatrix[0][1];
        cvaDvaRwr[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][5] = tempMatrix[0][0];
        cvaDvaRwr[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][7] = tempMatrix[0][1];
        cvaDvaRwr[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwr[i][9] = tempMatrix[0][0];
        cvaDvaRwr[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaRwr[i][11] = tempMatrix[0][1];
        cvaDvaRwr[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaRwr, "cvaDvaRwr.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Right Way Risk (G2++ Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;

    
    /*******************************************************************/
    /* RIGHT WAY RISK: G2++ QUASI-RANDOM                               */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaRwrSobol(exmpSwapTenors.size(),13);
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaRwrSobol[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        ext::shared_ptr<SobolRsg> sobolRsg(new SobolRsg(4*nTimeSteps, seed));
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrSobol[i][1] = tempMatrix[0][0];
        cvaDvaRwrSobol[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrSobol[i][3] = tempMatrix[0][1];
        cvaDvaRwrSobol[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrSobol[i][5] = tempMatrix[0][0];
        cvaDvaRwrSobol[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrSobol[i][7] = tempMatrix[0][1];
        cvaDvaRwrSobol[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, sobolRsg, nTimeSteps, nSamples,
                                        rho, g2Model, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrSobol[i][9] = tempMatrix[0][0];
        cvaDvaRwrSobol[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrSobol[i][11] = tempMatrix[0][1];
        cvaDvaRwrSobol[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaRwrSobol, "cvaDvaRwrSobol.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Right Way Risk (G2++ Quasi-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
    
    /*******************************************************************/
    /* RIGHT WAY RISK: HULL-WHITE PSEUDO-RANDOM                        */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    Matrix cvaDvaRwrHw(exmpSwapTenors.size(), 13);
    
    rhoHw[0][1]=rhoHw[1][0]=-0.99;
    
    for(int i=0;i<exmpSwapTenors.size();i++) {
        cvaDvaRwrHw[i][0] = years(exmpSwapTenors[i]);
        //one step per week
        Size nTimeSteps = 52*years(exmpSwapTenors[i]);
        
        tempSvapVec.clear();
        tempSvapVec.push_back(atmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrHw[i][1] = tempMatrix[0][0];
        cvaDvaRwrHw[i][2] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrHw[i][3] = tempMatrix[0][1];
        cvaDvaRwrHw[i][4] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(itmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrHw[i][5] = tempMatrix[0][0];
        cvaDvaRwrHw[i][6] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrHw[i][7] = tempMatrix[0][1];
        cvaDvaRwrHw[i][8] = qNorm95*tempMatrix[1][1];
        
        tempSvapVec.clear();
        tempSvapVec.push_back(otmSwapVec[i]);
        
        tempMatrix = mcCvaSwapPortfolio(tempSvapVec, uniformGenerator, nTimeSteps, nSamples,
                                        rhoHw, hwModel, cirCptModel, cirBankModel, cptRecovery, bankRecovery);
        
        cvaDvaRwrHw[i][9] = tempMatrix[0][0];
        cvaDvaRwrHw[i][10] = qNorm95*tempMatrix[1][0];
        cvaDvaRwrHw[i][11] = tempMatrix[0][1];
        cvaDvaRwrHw[i][12] = qNorm95*tempMatrix[1][1];
        
    }
    
    writeCsv(cvaDvaRwrHw, "cvaDvaRwrHw.csv");
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << "Right Way Risk (Hull-White Pseudo-Random): Run time = " << std::chrono::duration <double, std::ratio<60>> (diff).count() << " minutes" << std::endl;
}


#endif /* MonteCarlo_hpp */
