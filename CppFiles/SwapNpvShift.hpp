//
//  SwapNpvShift.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 31/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef SwapNpvShift_hpp
#define SwapNpvShift_hpp

#include "headers.hpp"


inline void swapNpvShift() {
    
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
    
    ext::shared_ptr<YieldTermStructure> yieldCurvePtr(new PiecewiseYieldCurve<Discount,Cubic>(anchorDate, rateHelperVec, floatingLegDayCounter));
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
    /* CREATING THE SWAP                                              */
    /*******************************************************************/
    
    Date startDate = calendar.advance(anchorDate,2, Days,ModifiedFollowing);
    Real nominal = 1000000.0;
    VanillaSwap::Type swapType = VanillaSwap::Payer;
    Rate tempFixedRate = 0.0;
    Spread swapSpread =0.0;
    
    
    Period swapTenor=30*Years;
    
    
    ext::shared_ptr<IborIndex> euriborIndex(new Euribor3M(yieldCurveHandle));
    
    euriborIndex->clearFixings(); //make sure that we don't have fixings from previous code
    
    
    VanillaSwap tempSwap=MakeVanillaSwap(swapTenor, euriborIndex,tempFixedRate)
    .withType(swapType)
    .withNominal(nominal)
    .withEffectiveDate(startDate)
    .withFloatingLegSpread(swapSpread)
    ;
    Real parRate = tempSwap.fairRate();
    
    VanillaSwap atmSwap=MakeVanillaSwap(swapTenor, euriborIndex,parRate)
    .withType(swapType)
    .withNominal(nominal)
    .withEffectiveDate(startDate)
    .withFloatingLegSpread(swapSpread)
    ;
    
    Matrix swapNpvShift(201,2);
    
    for(int i=0;i<201;i++) {
        Real x1=i*0.001-0.1;

    std::vector<Date> dates;
    std::vector<Real> discFactors;
    
    for(int j=0;j<300;j++) {
        
        dates.push_back(anchorDate+j*0.1*365);
        discFactors.push_back(g2Model->discountBond(0,j*0.1,{x1,0}));
    }
    

    ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> impliedTermStructure(
                                                                           new InterpolatedDiscountCurve<Cubic>(dates,discFactors,floatingLegDayCounter,calendar));
    
    impliedTermStructure->enableExtrapolation();
    
    
     Handle<YieldTermStructure> impliedHandle(impliedTermStructure);
    ext::shared_ptr<IborIndex> newIndex(new Euribor3M(impliedHandle));
    
    euriborIndex->clearFixings(); //make sure that we don't have fixings from previous code
    
    VanillaSwap newSwap=MakeVanillaSwap(swapTenor, newIndex,parRate)
    .withType(swapType)
    .withNominal(nominal)
    .withEffectiveDate(startDate)
    .withFloatingLegSpread(swapSpread)
    ;
    
        swapNpvShift[i][0]=x1;
        swapNpvShift[i][1]=newSwap.NPV();
    }
    writeCsv(swapNpvShift,"swapNpvShift.csv");
    
}

#endif /* SwapNpvShift_hpp */
