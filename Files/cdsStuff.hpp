//
//  cdsStuff.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 27/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef cdsStuff_hpp
#define cdsStuff_hpp

#include <stdio.h>

//#include <ql/instruments/makecds.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/experimental/credit/cdsoption.hpp>
#include <ql/experimental/credit/blackcdsoptionengine.hpp>

#include <ql/models/shortrate/calibrationhelpers/cdsoptionhelper.hpp>

using namespace QuantLib;

inline void cdsStuff() {
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
                                                                                            DateGeneration::CDS2015,Actual365Fixed(),
                                                                                            recovery_rate, yield_curve_handle)));
    }
    // Bootstrap hazard rates
    ext::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
    hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
                                                                            anchor_date, instruments, the_day_count));
    
    RelinkableHandle<DefaultProbabilityTermStructure> default_curve_handle;
    default_curve_handle.linkTo(hazardRateStructure);
    
    std::vector<std::pair<Date, Real> > hr_curve_data = hazardRateStructure->nodes();
    
    std::cout << "Calibrated hazard rate values: " << std::endl;
    for (auto& i : hr_curve_data) {
        std::cout << "hazard rate on " << i.first << " is " << i.second << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Some survival probability values: " << std::endl;
    std::cout << "1Y survival probability: "
    << io::percent(hazardRateStructure->survivalProbability(anchor_date +
                                                            1 * Years)) << std::endl;
    std::cout << "2Y survival probability: "
    << io::percent(hazardRateStructure->survivalProbability(anchor_date +
                                                            2 * Years)) << std::endl;
    
    std::cout << std::endl << std::endl;
    
    
    Period maturity=1*Years;
    Period length= 1*Years;
    Real recoveryRate = 0.4;
    Frequency paymentFrequency =Quarterly;
    Calendar calendar=TARGET();
    BusinessDayConvention paymentConvention =Following;
    BusinessDayConvention terminationDateConvention=Unadjusted;
    DateGeneration::Rule dateRule=DateGeneration::CDS2015;
    DayCounter dayCounter = Actual360();
    bool settlesAccrual = true;
    bool paysAtDefaultTime = true;
    Date protectionStart = Date();
    ext::shared_ptr< Claim >  claim   = ext::shared_ptr< Claim >();
    DayCounter lastPeriodDayCounter = DayCounter();
    bool rebatesAccrual = true;
    Date tradeDate = Date();
    Natural cashSettlementDays = 3;
    
    Handle<YieldTermStructure> discountCurve =yield_curve_handle;
    Handle<DefaultProbabilityTermStructure> defaultProbabilityCurve =default_curve_handle;
    
    
    Date startDate =the_calendar.advance(Settings::instance().evaluationDate(), maturity);
    Date endDate =the_calendar.advance(startDate,length);
    
    Schedule cds_schedule = MakeSchedule()
    .from(startDate)
    .to(endDate)
    .withFrequency(paymentFrequency)
    .withCalendar(calendar)
    .withConvention(paymentConvention)
    .withTerminationDateConvention(terminationDateConvention)
    .withRule(dateRule);
    
    Protection::Side side =Protection::Seller;
    
    //only running spread. 0.02 does not mean anything. We don't use this, but need to define it to find the fair spread.
    ext::shared_ptr<CreditDefaultSwap> temp(new CreditDefaultSwap(side, 1, 0.02, cds_schedule,paymentConvention,dayCounter,
                                                                  settlesAccrual,paysAtDefaultTime,protectionStart,claim,lastPeriodDayCounter, rebatesAccrual,tradeDate, cashSettlementDays));
    
    
    
    
    ext::shared_ptr<PricingEngine> cdsEngine(new MidPointCdsEngine(defaultProbabilityCurve, recoveryRate, discountCurve));
    
    temp->setPricingEngine(cdsEngine);
    
    
    std::cout << "CDS=" <<temp->NPV() << std::endl;
    
    std::cout << "Protection Start Date = " << temp->protectionStartDate() << std::endl;
    
    std::cout << "Fair Spread = " << temp->fairSpread() << std::endl;
    //only at the money
    Real forwardSpread =temp->fairSpread();
    
    ext::shared_ptr<CreditDefaultSwap> cds(new CreditDefaultSwap(side, 1, forwardSpread, cds_schedule,paymentConvention,dayCounter,
                                                                 settlesAccrual,paysAtDefaultTime,protectionStart,claim,lastPeriodDayCounter, rebatesAccrual,tradeDate, cashSettlementDays));
    
    cds->setPricingEngine(cdsEngine);
    std::cout << "CDS=" <<cds->NPV() << std::endl;
    
    ext::shared_ptr<EuropeanExercise> exercise(new EuropeanExercise(temp->protectionStartDate()));
    
    ext::shared_ptr<CdsOption> cdsOption;
    
    cdsOption=ext::make_shared<CdsOption>(cds, exercise);
    
    Handle<Quote> vol_handle=Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(0.15)));
    
    ext::shared_ptr<PricingEngine> black_pricing_engine(new BlackCdsOptionEngine(default_curve_handle, recovery_rate, yield_curve_handle,vol_handle));
    
    cdsOption->setPricingEngine(black_pricing_engine);
    
    std::cout << "CDS Option=" <<cdsOption->NPV() << std::endl;
    
    CdsOptionHelper test(1*Years,1*Years, vol_handle, 0.4, default_curve_handle, yield_curve_handle);
    test.setPricingEngine(black_pricing_engine);
    
    std::cout << test.modelValue() <<std::endl;
    
    
}

#endif /* cdsStuff_hpp */
