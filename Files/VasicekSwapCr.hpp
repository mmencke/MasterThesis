//
//  VasicekSwapCr.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 24/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef VasicekSwapCr_hpp
#define VasicekSwapCr_hpp

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

#include <ql/pricingengines/swaption/jamshidianswaptionengine.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>
#include <ql/termstructures/credit/interpolatedhazardratecurve.hpp>
#include <ql/pricingengines/swap/cvaswapengine.hpp>

using namespace QuantLib;

inline void vasicekSwapCr() {
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
    Period Swap_Length = 30 *Years;
    
    VanillaSwap myNewSwap = MakeVanillaSwap(Swap_Length, euriborIndex,Fixed_Rate)
    .withType(Swap_Type)
    .withNominal(Nominal)
    .withEffectiveDate(Start_Date)
    .withFloatingLegSpread(Swap_Spread)
    ;
    
    
    Schedule newFloat = myNewSwap.floatingSchedule();
    
    
    Schedule newFixed = myNewSwap.fixedSchedule();
    
    Real kappa = 0.5;
    Real theta =0.03;
    Real sigma = 0.1; //absolute volatility
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
    
    ext::shared_ptr<InterpolatedZeroCurve<BackwardFlat>> the_curve_ptr(new InterpolatedZeroCurve<BackwardFlat>(the_dates,the_rates, Float_Day_Count_Basis, This_Calendar));
    
    the_handle.linkTo(the_curve_ptr);
    
    Real Anchor_Fixing = euriborIndex->forecastFixing(Anchor_Date);//fixing to forecast so consistent with zero curve
    
    euriborIndex->addFixing(Anchor_Date, Anchor_Fixing);
    
    
    ext::shared_ptr<Vasicek> vasicek_model_ptr(new Vasicek(vasicek_model)); //clone the vasicek model to ptr
    
    ext::shared_ptr<JamshidianSwaptionEngine> swaption_engine_ptr(new JamshidianSwaptionEngine(vasicek_model_ptr, the_handle));
    //ext::shared_ptr<BlackSwaptionEngine> swaption_engine_ptr(new BlackSwaptionEngine(the_handle,0.15));
    
    Size defaultTenors[] = {0, 12, 36, 60, 84, 120, 180, 240, 300,
        360};// months
    // Three risk levels:
    Real intensitiesLow[] = {0.0534, 0.0534, 0.0564, 0.06, 0.0614, 0.0696,
        0.0696, 0.0696, 0.0696, 0.0696, 0.0696};
    // Recovery rates:
    Real ctptyRRLow = 0.3;
    
    std::vector<Date> defaultTSDates;
    std::vector<Real> intesitiesVLow, intesitiesVMedium, intesitiesVHigh;
    
    for(Size i=0; i<sizeof(defaultTenors)/sizeof(Size); i++) {
        defaultTSDates.push_back(TARGET().advance(Anchor_Date,
                                                  Period(defaultTenors[i], Months)));
        intesitiesVLow.push_back(intensitiesLow[i]);
    }
    std::vector<Handle<DefaultProbabilityTermStructure> >
    defaultIntensityTS;
    
    defaultIntensityTS.push_back(Handle<DefaultProbabilityTermStructure>(
                                                                         ext::shared_ptr<DefaultProbabilityTermStructure>(
                                                                                                                          new InterpolatedHazardRateCurve<BackwardFlat>(
                                                                                                                                                                        defaultTSDates,
                                                                                                                                                                        intesitiesVLow,
                                                                                                                                                                        Actual360(),
                                                                                                                                                                        TARGET()))));
    defaultIntensityTS[0]->enableExtrapolation();
    Handle<PricingEngine> the_pricing_handle(swaption_engine_ptr);
    
    ext::shared_ptr<PricingEngine> ctptySwapCvaLow =
    ext::make_shared<CounterpartyAdjSwapEngine>(
                                                the_handle,
                                                the_pricing_handle,
                                                defaultIntensityTS[0],
                                                ctptyRRLow
                                                );
    
    ext::shared_ptr<VanillaSwap> riskySwapPtr(new VanillaSwap(myNewSwap)); //clone the swap to ptr
    
    riskySwapPtr->setPricingEngine(ctptySwapCvaLow);
    
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the swap without credit risk is" % " : " % myNewSwap.NPV() << std::endl;
    std::cout << boost::format("%-30s %-3s%-9s") % "The par swap rate without credit risk is" % " : " % myNewSwap.fairRate() << std::endl;
    
    std::cout << boost::format("%-30s %-3s %-9s") % "The NPV of the swap with credit risk is" % " : " % riskySwapPtr->NPV() << std::endl;
    std::cout << boost::format("%-30s %-3s%-9s") % "The par swap rate with credit risk is" % " : " % riskySwapPtr->fairRate() << std::endl;
    
    Real bp_adjustment =(myNewSwap.fairRate()-riskySwapPtr->fairRate())*10000;
    
    std::cout << boost::format("%-30s %-3s%-9s") % "BP Adjustment = " % " : " % bp_adjustment << std::endl;
}

#endif /* VasicekSwapCr_hpp */
