//
//  McCirCds.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 27/04/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef McCirCds_hpp
#define McCirCds_hpp

#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>

#include <ql/math/distributions/chisquaredistribution.hpp>

#include <ql/math/randomnumbers/mt19937uniformrng.hpp>

#include <ql/experimental/credit/mccircdsoptionengine.hpp>

#include <ql/methods/montecarlo/montecarlomodel.hpp>
#include <ql/models/shortrate/onefactormodels/extendedcoxingersollross.hpp>

using namespace QuantLib;

inline void mcCirCds() {
    
    
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
    hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(anchor_date, instruments, the_day_count));
    
    RelinkableHandle<DefaultProbabilityTermStructure> default_curve_handle;
    default_curve_handle.linkTo(hazardRateStructure);

    
    Period maturity=1*Years;
    Period length= 10*Years;
    Real recoveryRate = 0.4;
    Frequency paymentFrequency =Quarterly;
    Calendar calendar=TARGET();
    BusinessDayConvention paymentConvention =Following;
    BusinessDayConvention terminationDateConvention=Unadjusted;
    DateGeneration::Rule dateRule=DateGeneration::CDS2015;
    DayCounter dayCounter = Actual360();
    
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
    ext::shared_ptr<CreditDefaultSwap> cds(new CreditDefaultSwap(side, 1, 0.02, cds_schedule,paymentConvention,dayCounter));
    
    
    
    
    ext::shared_ptr<PricingEngine> cdsEngine(new MidPointCdsEngine(defaultProbabilityCurve, recoveryRate, discountCurve));
    
    cds->setPricingEngine(cdsEngine);
    
    for(auto x : cds_schedule.dates()) {
        std::cout << x << std::endl;
    }
    
    //Settings::instance().evaluationDate() =Date(22,March,2010);
    
    startDate = Settings::instance().evaluationDate();
    Date exerciseDate = cds->protectionStartDate();
    Date maturityDate = cds->protectionEndDate();
    
    
    BigNatural seed = 1;
    Size n_time_steps = 1; //only one time step
    Size n_samples = 100000;
    Statistics statistics_accumulator;
    
    bool brownian_bridge = false;
    bool antithetic_variate = false;
    
    //from page 125 in Brigo, Morini and Pallavicini
    Real kappa=0.4;
    Real theta = 0.026;
    Real sigma = 0.14;
    Real lambda0 = 0.0165;
    
    //This is only used for determinng how long to simulate
    Real t=Actual360().yearFraction(startDate,exerciseDate);
    Real T=Actual360().yearFraction(startDate,maturityDate);
    /*
     Real c=(4*kappa)/(sigma*sigma*(1-std::exp(-kappa*t)));
     Real nu=(4*kappa*theta)/(sigma*sigma);
     Real eta=c*lambda0*std::exp(-kappa*t);
     
     MersenneTwisterUniformRng rng(seed);
     
     InverseNonCentralCumulativeChiSquareDistribution chi2(nu, eta,100);
     
     std::cout << "Random Chi2=" << chi2(rng.nextReal())/c << std::endl;
     */
    TimeGrid time_grid(t,n_time_steps);
    
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(n_time_steps, seed);
    
    ext::shared_ptr<CoxIngersollRoss> cir_model(new CoxIngersollRoss(lambda0, theta, kappa, sigma));
    
    
    //ext::shared_ptr<StochasticProcess1D> cir_process(new CoxIngersollRossProcess(kappa, sigma, lambda0, theta,CoxIngersollRossProcess::Exact));
    
    ext::shared_ptr<StochasticProcess1D> cir_process = cir_model->dynamics()->process();
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cir_rsg_ptr(new SingleVariate<PseudoRandom>::path_generator_type(cir_process,time_grid, rsg,false));
    
    Matrix cir_single_step(n_samples,1);
    
    
    auto start = std::chrono::steady_clock::now();
    for(int i=0;i<n_samples;i++) {
        Path path = (cir_rsg_ptr->next()).value;
        cir_single_step[i][0] =path[n_time_steps];
        //cir_single_step[i][0] =chi2(rng.nextReal())/c;
    }
    auto end = std::chrono::steady_clock::now();
    
    auto diff = end - start;
    
    std::cout << "run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    
    writeCsv(cir_single_step, "cir_single_step.csv");
    
    
    
    Real n=(T-t)*4; //we want 4 timesteps per year (once per quarter)
    Real dt=(T-t)/n;
    
    std::cout << "protection start" << exerciseDate << std::endl;
    
    std::cout << "protection end" << maturityDate << std::endl;
    
    std::cout << "n=" << n << std::endl;
    std::cout << "dt=" << dt << std::endl;
    std::vector<Date> dates;
    std::vector<Probability> survivalProbabilities;
    
    dates.push_back(exerciseDate);
    survivalProbabilities.push_back(1);
    
    for(int i=1;i<=n;i++) {
        dates.push_back(exerciseDate+i*dt*365);
        survivalProbabilities.push_back(cir_model->discountBond(t,t+i*dt,cir_single_step[0][0]));
        
    }
    
    ext::shared_ptr<InterpolatedSurvivalProbabilityCurve<Cubic>> x(new InterpolatedSurvivalProbabilityCurve<Cubic>(dates, survivalProbabilities, Actual360()));
    x->enableExtrapolation(true);
    int m= 100;
    
    Matrix test(m,2);
    for(int i=0;i<m;i++) {
        test[i][0]=i*0.1;
        test[i][1]=x->survivalProbability(test[i][0],true);
    }
    writeCsv(test,"test.csv");
    
    //Settings::instance().evaluationDate() = exerciseDate;
    ext::shared_ptr<PricingEngine> cdsEngine2(new MidPointCdsEngine(Handle<DefaultProbabilityTermStructure>(x), recoveryRate, discountCurve));
    
    cds->setPricingEngine(cdsEngine2);
    
    std::cout << cds->NPV() << std::endl;
    
    /*ext::shared_ptr<PricingEngine> cdsEngine3(new McCirCdsEngine(Handle<DefaultProbabilityTermStructure>(x), recoveryRate, discountCurve));
    
    cds->setPricingEngine(cdsEngine3);
    
    std::cout << cds->NPV() << std::endl;*/
    
    
    ext::shared_ptr<EuropeanExercise> exercise(new EuropeanExercise(cds->protectionStartDate()));
    
    ext::shared_ptr<CdsOption> cdsOption;
    
    cdsOption=ext::make_shared<CdsOption>(cds, exercise);
    
    
    ext::shared_ptr<PricingEngine> cdsEngine3(new McCirCdsOptionEngine(cir_model, 100, seed, recoveryRate, discountCurve));
    cdsOption->setPricingEngine(cdsEngine3);
    std::cout << cdsOption->NPV() << std::endl;
    
    
    ExtendedCoxIngersollRoss xCirModel(yield_curve_handle);
    
    
    for(auto x : xCirModel.params()) {
        std::cout << x << std::endl;
    }
    
    
    
}


#endif /* McCirCds_hpp */
