//
//  SimpleSwapPricing.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef SimpleSwapPricing_hpp
#define SimpleSwapPricing_hpp

#include <iostream>

#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>

#include <ql/indexes/ibor/euribor.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/instruments/makevanillaswap.hpp>

#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>

#include <boost/format.hpp>

using namespace QuantLib;

inline void simpleSwapPricing() {
    std::cout << "Here we are pricing a plain vanilla 3m/1y swap." << std::endl << std::endl;
    
    Date Anchor_Date(16,February,2021);
    
    TARGET This_Calendar;
    
    Date Start_Date = This_Calendar.advance(Anchor_Date,2, Days,ModifiedFollowing);
    
    Actual360 Float_Day_Count_Basis;
    
    Real Nominal = 1000000.0;
    
    VanillaSwap::Type Swap_Type = VanillaSwap::Payer;
    
    Rate Fixed_Rate = 0.05;

    RelinkableHandle<YieldTermStructure> myHandle;
    
    ext::shared_ptr<IborIndex> euriborIndex(new Euribor3M(myHandle));
    
    Spread Swap_Spread =0.0;
    
    
    VanillaSwap myNewSwap = MakeVanillaSwap(1*Years, euriborIndex,Fixed_Rate)
    .withType(Swap_Type)
    .withNominal(Nominal)
    .withEffectiveDate(Start_Date)
    .withFloatingLegSpread(Swap_Spread)
    ;
    
    Schedule newFloat = myNewSwap.floatingSchedule();
    
    std::cout << "Number of floating payments: " << newFloat.size() << std::endl;
    std::cout << "Dates are ouput below:" << std::endl;
    
    for(auto x : newFloat.dates()) {
        std::cout << "   " << x << std::endl;
    }
    
    Schedule newFixed = myNewSwap.fixedSchedule();
    
    std::cout << std::endl << "Number of fixed payments: " << newFixed.size() << std::endl;
    std::cout << "Dates are ouput below:" << std::endl;
    for(auto x : newFixed.dates()) {
        std::cout << "   " << x << std::endl;
    }
    
    std::vector<Date> the_dates = newFloat.dates();
    std::vector<Rate> the_rates = {0.01,0.02,0.03,0.04,0.05};
    
    ext::shared_ptr<InterpolatedZeroCurve<BackwardFlat>> the_curve_ptr(new InterpolatedZeroCurve<BackwardFlat>(the_dates,the_rates, Float_Day_Count_Basis, This_Calendar));
    
    myHandle.linkTo(the_curve_ptr);
    
    
    Real Anchor_Fixing = euriborIndex->forecastFixing(Anchor_Date);//fixing to forecast so consistent with zero curve
    
    euriborIndex->addFixing(Anchor_Date, Anchor_Fixing);
    
    
    std::cout << std::endl << "We have the valuation on " << myNewSwap.valuationDate() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the swap is" % " : " % myNewSwap.NPV() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the fixed leg is" % " : " % myNewSwap.legNPV(0) << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the floating leg is" % " : " % myNewSwap.legNPV(1) << std::endl;
    std::cout << boost::format("%-30s %-3s%-9s") % "The par swap rate is" % " : " % myNewSwap.fairRate() << std::endl;
    
    
}
#endif /* SimpleSwapPricing_hpp */

