//
//  VasicekCleanSwaphpp
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef VasicekCleanSwap_hpp
#define VasicekCleanSwap_hpp

#include <iostream>
#include <boost/format.hpp>

#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>

#include <ql/indexes/ibor/euribor.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/instruments/makevanillaswap.hpp>

#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>

#include <ql/models/shortrate/onefactormodels/vasicek.hpp>

#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>

#include "utility.hpp"

using namespace QuantLib;

inline void vasicekSwapNoCr() {
    std::cout << "Here we are pricing the same swap as before in Vasicek." << std::endl << std::endl;
    
    Date Anchor_Date(16,February,2021);
    Settings::instance().evaluationDate() = Anchor_Date;
    
    TARGET This_Calendar;
    
    Date Start_Date = This_Calendar.advance(Anchor_Date,2, Days,ModifiedFollowing);
    
    Actual360 Float_Day_Count_Basis;
    
    Real Nominal = 1000000.0;
    
    VanillaSwap::Type Swap_Type = VanillaSwap::Receiver;
    
    Rate Fixed_Rate = 0.05;

    RelinkableHandle<YieldTermStructure> the_handle;
    
    ext::shared_ptr<IborIndex> euriborIndex(new Euribor3M(the_handle));
    
    euriborIndex->clearFixings(); //make sure that we don't have fixings from previous code
    
    Spread Swap_Spread =0.0;
    
    
    VanillaSwap myNewSwap = MakeVanillaSwap(1*Years, euriborIndex,Fixed_Rate)
    .withType(Swap_Type)
    .withNominal(Nominal)
    .withEffectiveDate(Start_Date)
    .withFloatingLegSpread(Swap_Spread)
    ;
    
    
    Schedule newFloat = myNewSwap.floatingSchedule();
    
    
    Schedule newFixed = myNewSwap.fixedSchedule();
    
    Real kappa = 0.5;
    Real theta =0.03;
    Real sigma = 0.01; //absolute volatility
    Real r0 = 0.01;
    
    Vasicek vasicek_model(r0, kappa, theta,sigma);
    
    std::vector<Date> the_dates = newFloat.dates();
    
    Size n =the_dates.size();
    std::vector<Rate> the_rates;
    
    std::vector<Time> year_fractions;
    
    for(int i=0;i<n;i++) {
        year_fractions.push_back(Float_Day_Count_Basis.yearFraction(Anchor_Date,the_dates[i]));
    }
    
    for(int i=0;i<n;i++) {
        the_rates.push_back(-log(vasicek_model.discount(year_fractions[i]))/year_fractions[i]);
    }
    
    std::cout << "We have rates for dates" << std::endl;
    
    for(int i=0;i<n;i++) {
        std::cout << boost::format("%-3s %-10s %-3s %-9s") % "   " % the_dates[i] % " : " % the_rates[i];
        std::cout << std::endl;
    }
    
    //the_rates = {0.01,0.02,0.03,0.04,0.05};
    
    ext::shared_ptr<InterpolatedZeroCurve<BackwardFlat>> the_curve_ptr(new InterpolatedZeroCurve<BackwardFlat>(the_dates,the_rates, Float_Day_Count_Basis, This_Calendar));
    
    the_handle.linkTo(the_curve_ptr);
    
    Real Anchor_Fixing = euriborIndex->forecastFixing(Anchor_Date);//fixing to forecast so consistent with zero curve
    
    euriborIndex->addFixing(Anchor_Date, Anchor_Fixing);
    
    
    
    std::cout << std::endl << "We have the valuation on " << myNewSwap.valuationDate() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the swap is" % " : " % myNewSwap.NPV() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the fixed leg is" % " : " % myNewSwap.legNPV(0) << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the floating leg is" % " : " % myNewSwap.legNPV(1) << std::endl;
    std::cout << boost::format("%-30s %-3s%-9s") % "The par swap rate is" % " : " % myNewSwap.fairRate() << std::endl;

    std::cout << std::endl;
    
    Date New_Anchor = Date(15,November,2021);
    Settings::instance().evaluationDate() = New_Anchor;
    
    euriborIndex->addFixing(Date(16,August,2021), 0.0087344);
    
    int j=0;
    std::vector<Date> the_new_dates;
    
    the_new_dates.push_back(New_Anchor);
    
    for(int i =0; i<n;i++) {
        if(the_dates[i]>New_Anchor) {
            the_new_dates.push_back(the_dates[i]);
            j+=1;
        }
    }
    
    n = the_new_dates.size()-1;
    
    std::vector<Rate> the_new_rates;
    
    the_new_rates.push_back(vasicek_model.r0());
    
    std::vector<Time> new_year_fractions;
    
    for(int i=0;i<n;i++) {
        new_year_fractions.push_back(Float_Day_Count_Basis.yearFraction(New_Anchor,the_new_dates[i+1]));
    }
    
    for(int i=0;i<n;i++) {
        the_new_rates.push_back(-log(vasicek_model.discount(new_year_fractions[i]))/new_year_fractions[i]);
    }
    
    std::cout << "We have rates for dates" << std::endl;
    
    for(int i=0;i<n+1;i++) {
        std::cout << boost::format("%-3s %-10s %-3s %-9s") % "   " % the_new_dates[i] % " : " % the_new_rates[i];
        std::cout << std::endl;
    }
    
    ext::shared_ptr<InterpolatedZeroCurve<BackwardFlat>> the_new_curve_ptr(new InterpolatedZeroCurve<BackwardFlat>(the_new_dates,the_new_rates, Float_Day_Count_Basis, This_Calendar));
    
    the_handle.linkTo(the_new_curve_ptr);
    
    std::cout << std::endl << "We have the valuation on " << myNewSwap.valuationDate() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the swap is" % " : " % myNewSwap.NPV() << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the fixed leg is" % " : " % myNewSwap.legNPV(0) << std::endl;
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the floating leg is" % " : " % myNewSwap.legNPV(1) << std::endl;
    std::cout << boost::format("%-30s %-3s%-9s") % "The par swap rate is" % " : " % myNewSwap.fairRate() << std::endl;

    std::cout << std::endl << "some test:"<< std::endl;
    
    
    Real dt = 0.35;
    
    std::cout << Anchor_Date + dt*365 << std::endl;
    
    std::cout << "Expiry date=" <<the_new_dates[1] << std::endl;
    
    
    
    
    ext::shared_ptr<VanillaSwap> mySwapPtr(new VanillaSwap(myNewSwap)); //clone the swap to ptr

    
    Real blackVol=0.15;
    Swaption mySwaption(mySwapPtr,ext::make_shared<EuropeanExercise>(the_new_dates[1]));
    mySwaption.setPricingEngine(ext::make_shared<BlackSwaptionEngine>(the_handle, blackVol));
    
    
    //std::cout << "The Swaption NPV is " << mySwaption.NPV() << std::endl;
}
#endif /* VasicekCleanSwap_hpp */

