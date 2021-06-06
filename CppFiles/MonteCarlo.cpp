//
//  McCva.cpp
//  Master Thesis
//
//  Created by Magnus Mencke on 11/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#include "MonteCarlo.hpp"


Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap>> swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<G2> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery) {
    
    Date anchorDate = Settings::instance().evaluationDate();

    Date longestMaturity =swapPortfolio[0]->maturityDate();
    for(int i=0;i<swapPortfolio.size();i++) {
        if(longestMaturity<swapPortfolio[i]->maturityDate()) {
            longestMaturity=swapPortfolio[i]->maturityDate();
        }
    }
    
    Actual365Fixed dummyDayCounter;
    TARGET dummyCalendar;
    
    TimeGrid timeGrid(dummyDayCounter.yearFraction(anchorDate,longestMaturity),nTimeSteps);
    
    InverseCumulativeRng<MersenneTwisterUniformRng,InverseCumulativeNormal> rng(*uniformGenerator);
    
    Matrix rhoChol = CholeskyDecomposition(rho);
    
    Natural nNormals = 4;//we need to have 4 normals for every step
    
    Real pathWeight=1;
    
    Statistics cvaAccumulator;
    Statistics dvaAccumulator;
    
    ext::shared_ptr<StochasticProcess> interestRateProcess = interestRate->dynamics()->process();
    ext::shared_ptr<StochasticProcess> ctptyIntensityProcess = ctptyIntensity->dynamics()->process();
    ext::shared_ptr<StochasticProcess> invstIntensityProcess = invstIntensity->dynamics()->process();
    
    Matrix cvaDvaSample(nSamples,4);
    
    for(int l=0;l<nSamples;l++) {
        Real ctptyU =1-uniformGenerator->nextReal();
        Real invstU =1-uniformGenerator->nextReal();
        
        Array x=interestRateProcess->initialValues();
        Array y=ctptyIntensityProcess->initialValues();
        Array z=invstIntensityProcess->initialValues();
        
        
        
        Real bankAccount=1;
        Real ctptySurvProb=1;
        Real invstSurvProb=1;
        
        bool ctptyDefault=false;
        bool invstDefault=false;
        
        Real portfolioNpv=0;
        Time tau;
        
        Real cva,dva;
        
        for(int i=0;i<nTimeSteps-1;i++) {//default on the last timestep does not count as a default
            Time t = timeGrid[i];
            Time dt =timeGrid.dt(i);
            Array dw(nNormals);
            for(int j=0;j<nNormals;j++) {
                dw[j]=rng.next().value;
            }
            Array dwCor(nNormals);
            dwCor=rhoChol*dw;
            
            Array tempX={dwCor[0],dwCor[1]};
            Array tempY={dwCor[2]};
            Array tempZ={dwCor[3]};
            
            x=interestRateProcess->evolve(t,x,dt,tempX);
            y=ctptyIntensityProcess->evolve(t,y,dt,tempY);
            z=invstIntensityProcess->evolve(t,z,dt,tempZ);
            
            bankAccount/=interestRate->discountBond(t,t+dt,x);
            ctptySurvProb*=ctptyIntensity->discountBond(t,t+dt,y);
            invstSurvProb*=invstIntensity->discountBond(t,t+dt,z);
            
            if(ctptyU>=ctptySurvProb) {
                ctptyDefault=true;
                tau=t;
                break;
            }
            
            if(invstU>=invstSurvProb){
                invstDefault=true;
                tau=t;
                break;
            }
            
        }
        if(ctptyDefault||invstDefault) {
            Date defaultDate=anchorDate+tau*365;
            Settings::instance().evaluationDate()=defaultDate;
            
            std::vector<Date> dates;
            std::vector<Real> discFactors;
            
            for(int j=timeGrid.index(tau);j<nTimeSteps;j++) {
                dates.push_back(anchorDate+timeGrid[j]*365);
                discFactors.push_back(interestRate->discountBond(tau,timeGrid[j],x));
            }
            
            ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> impliedTermStructure(
                                                                                   new InterpolatedDiscountCurve<Cubic>(dates,discFactors,dummyDayCounter,dummyCalendar));
            
            impliedTermStructure->enableExtrapolation();
            
            Handle<YieldTermStructure> impliedHandle(impliedTermStructure);
            for(int k=0;k<swapPortfolio.size();k++) {
                
                ext::shared_ptr<IborIndex> newIndex(
                                                    new IborIndex(swapPortfolio[k]->iborIndex()->familyName(),
                                                                  swapPortfolio[k]->iborIndex()->tenor(),
                                                                  swapPortfolio[k]->iborIndex()->fixingDays(),
                                                                  swapPortfolio[k]->iborIndex()->currency(),
                                                                  swapPortfolio[k]->iborIndex()->fixingCalendar(),
                                                                  swapPortfolio[k]->iborIndex()->businessDayConvention(),
                                                                  swapPortfolio[k]->iborIndex()->endOfMonth(),
                                                                  swapPortfolio[k]->iborIndex()->dayCounter(),
                                                                  impliedHandle));
                newIndex->clearFixings();//to remove fixings from other paths
                VanillaSwap temp(swapPortfolio[k]->type(),
                                 swapPortfolio[k]->nominal(),
                                 swapPortfolio[k]->fixedSchedule(),
                                 swapPortfolio[k]->fixedRate(),
                                 swapPortfolio[k]->fixedDayCount(),
                                 swapPortfolio[k]->floatingSchedule(),
                                 newIndex,
                                 swapPortfolio[k]->spread(),
                                 swapPortfolio[k]->floatingDayCount(),
                                 swapPortfolio[k]->paymentConvention());
                //discounting using the forward curve
                Handle<YieldTermStructure> disc =newIndex->forwardingTermStructure();
                ext::shared_ptr<PricingEngine> engine(new
                                                      DiscountingSwapEngine(disc));
                
                temp.setPricingEngine(engine);
                
                Date dateBefore =temp.floatingSchedule().previousDate(defaultDate);
                Date dateAfter1 =temp.floatingSchedule().nextDate(defaultDate);
                Date dateAfter2 =temp.floatingSchedule().nextDate(dateAfter1);
                
                //approximate the last fixing as the current simple forward rate
                Rate forwardRate1 = impliedTermStructure->forwardRate(defaultDate,dateAfter1, temp.floatingDayCount(), Simple);
                Rate forwardRate2 = impliedTermStructure->forwardRate(dateAfter1,dateAfter2, temp.floatingDayCount(), Simple);
                
                if(dateBefore!=Date()) {//if not null date
                    Date fixing1 =newIndex->fixingDate(dateBefore);
                    newIndex->addFixing(fixing1,forwardRate1);
                }
                
                Date fixing2 =newIndex->fixingDate(dateAfter1);
                if(fixing2<defaultDate) {
                    newIndex->addFixing(fixing2,forwardRate2);
                }
                
                portfolioNpv+=temp.NPV();
            }
            
            
            if(ctptyDefault) {
                if(portfolioNpv>0) {
                    cva=(1-ctptyRecovery)*portfolioNpv/bankAccount;
                } else {
                    cva=0;
                }
            }
            
            if(invstDefault) {
                if(portfolioNpv<0) {
                    dva=-(1-invstRecovery)*portfolioNpv/bankAccount;
                } else {
                    dva=0;
                }
            }
            
            
        }
        
        if (!ctptyDefault) {
            cva=0;
        }
        if (!invstDefault) {
            dva=0;
        }
        
        cvaAccumulator.add(cva,pathWeight);
        dvaAccumulator.add(dva,pathWeight);
        
        cvaDvaSample[l][0]=cva;
        cvaDvaSample[l][1]=dva;
        cvaDvaSample[l][2]=x[0];
        cvaDvaSample[l][3]=x[1];
        
        
    }
    
    writeCsv(cvaDvaSample,"cvaDvaSample.csv");
    
    //reset evaluationdate
    Settings::instance().evaluationDate()=anchorDate;
    
    Matrix results(2,2);
    results[0][0]=cvaAccumulator.mean();
    results[1][0]=cvaAccumulator.errorEstimate();
    results[0][1]=dvaAccumulator.mean();
    results[1][1]=dvaAccumulator.errorEstimate();
    
    return results;
}




Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap>> swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, ext::shared_ptr<SobolRsg> sobolRsg, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<G2> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery) {
    
    Date anchorDate = Settings::instance().evaluationDate();
    
    Date longestMaturity =swapPortfolio[0]->maturityDate();
    for(int i=0;i<swapPortfolio.size();i++) {
        if(longestMaturity<swapPortfolio[i]->maturityDate()) {
            longestMaturity=swapPortfolio[i]->maturityDate();
        }
    }
    
    Actual365Fixed dummyDayCounter;
    TARGET dummyCalendar;
    
    TimeGrid timeGrid(dummyDayCounter.yearFraction(anchorDate,longestMaturity),nTimeSteps);
    
    
    
    std::vector<ext::shared_ptr<StochasticProcess1D>> processVec;
    
    processVec.push_back(interestRate->dynamics()->xProcess());
    processVec.push_back(interestRate->dynamics()->yProcess());
    processVec.push_back(ctptyIntensity->dynamics()->process());
    processVec.push_back(invstIntensity->dynamics()->process());
    
    
    Matrix rhoAdj=rho;
    
    rhoAdj[0][1]=rhoAdj[1][0]=interestRate->rho();
    
    ext::shared_ptr<StochasticProcess> processArray(new StochasticProcessArray(processVec,rhoAdj));

    
    
    InverseCumulativeNormal quantileNormal;
    LowDiscrepancy::rsg_type rsg =InverseCumulativeRsg<QuantLib::SobolRsg, QuantLib::InverseCumulativeNormal>(*sobolRsg,quantileNormal);
    
    bool brownianBridge=false;
    ext::shared_ptr<MultiVariate<LowDiscrepancy>::path_generator_type> multiPathGenerator(new MultiVariate<LowDiscrepancy>::path_generator_type(processArray,timeGrid,rsg, brownianBridge));
    
    Real pathWeight=1;
    
    Statistics cvaAccumulator;
    Statistics dvaAccumulator;
    
    Matrix cvaDvaSample(nSamples,4);
    
    for(int l=0;l<nSamples;l++) {
        Real ctptyU =1-uniformGenerator->nextReal();
        Real invstU =1-uniformGenerator->nextReal();
        
        
        MultiPath thePath = multiPathGenerator->next().value;
        
        Array x={thePath[0][0],thePath[1][0]};
        Array y={thePath[2][0]};
        Array z={thePath[3][0]};
        
        
        Real bankAccount=1;
        Real ctptySurvProb=1;
        Real invstSurvProb=1;
        
        bool ctptyDefault=false;
        bool invstDefault=false;
        
        Real portfolioNpv=0;
        Time tau;
        
        Real cva,dva;
        
        for(int i=0;i<nTimeSteps-1;i++) {//default on the last timestep does not count as a default
            
            Time t = timeGrid[i];
            Time dt =timeGrid.dt(i);
            
            
            x={thePath[0][i+1],thePath[1][i+1]};
            y={thePath[2][i+1]};
            z={thePath[3][i+1]};
            
            bankAccount/=interestRate->discountBond(t,t+dt,x);
            ctptySurvProb*=ctptyIntensity->discountBond(t,t+dt,y);
            invstSurvProb*=invstIntensity->discountBond(t,t+dt,z);
            
            if(ctptyU>=ctptySurvProb) {
                ctptyDefault=true;
                tau=t;
                break;
            }
            
            if(invstU>=invstSurvProb){
                invstDefault=true;
                tau=t;
                break;
            }
            
        }
        if(ctptyDefault||invstDefault) {
            Date defaultDate=anchorDate+tau*365;
            Settings::instance().evaluationDate()=defaultDate;
            
            std::vector<Date> dates;
            std::vector<Real> discFactors;
            
            for(int j=timeGrid.index(tau);j<nTimeSteps;j++) {
                dates.push_back(anchorDate+timeGrid[j]*365);
                discFactors.push_back(interestRate->discountBond(tau,timeGrid[j],x));
            }
            
            ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> impliedTermStructure(
                                                                                   new InterpolatedDiscountCurve<Cubic>(dates,discFactors,dummyDayCounter,dummyCalendar));
            
            impliedTermStructure->enableExtrapolation();
            
            Handle<YieldTermStructure> impliedHandle(impliedTermStructure);
            for(int k=0;k<swapPortfolio.size();k++) {
                
                ext::shared_ptr<IborIndex> newIndex(
                                                    new IborIndex(swapPortfolio[k]->iborIndex()->familyName(),
                                                                  swapPortfolio[k]->iborIndex()->tenor(),
                                                                  swapPortfolio[k]->iborIndex()->fixingDays(),
                                                                  swapPortfolio[k]->iborIndex()->currency(),
                                                                  swapPortfolio[k]->iborIndex()->fixingCalendar(),
                                                                  swapPortfolio[k]->iborIndex()->businessDayConvention(),
                                                                  swapPortfolio[k]->iborIndex()->endOfMonth(),
                                                                  swapPortfolio[k]->iborIndex()->dayCounter(),
                                                                  impliedHandle));
                newIndex->clearFixings();//to remove fixings from other paths
                VanillaSwap temp(swapPortfolio[k]->type(),
                                 swapPortfolio[k]->nominal(),
                                 swapPortfolio[k]->fixedSchedule(),
                                 swapPortfolio[k]->fixedRate(),
                                 swapPortfolio[k]->fixedDayCount(),
                                 swapPortfolio[k]->floatingSchedule(),
                                 newIndex,
                                 swapPortfolio[k]->spread(),
                                 swapPortfolio[k]->floatingDayCount(),
                                 swapPortfolio[k]->paymentConvention());
                //discounting using the forward curve
                Handle<YieldTermStructure> disc =newIndex->forwardingTermStructure();
                ext::shared_ptr<PricingEngine> engine(new
                                                      DiscountingSwapEngine(disc));
                
                temp.setPricingEngine(engine);
                
                Date dateBefore =temp.floatingSchedule().previousDate(defaultDate);
                Date dateAfter1 =temp.floatingSchedule().nextDate(defaultDate);
                Date dateAfter2 =temp.floatingSchedule().nextDate(dateAfter1);
                
                //approximate the last fixing as the current simple forward rate
                Rate forwardRate1 = impliedTermStructure->forwardRate(defaultDate,dateAfter1, temp.floatingDayCount(), Simple);
                Rate forwardRate2 = impliedTermStructure->forwardRate(dateAfter1,dateAfter2, temp.floatingDayCount(), Simple);
                
                if(dateBefore!=Date()) {//if not null date
                    Date fixing1 =newIndex->fixingDate(dateBefore);
                    newIndex->addFixing(fixing1,forwardRate1);
                }
                
                Date fixing2 =newIndex->fixingDate(dateAfter1);
                if(fixing2<defaultDate) {
                    newIndex->addFixing(fixing2,forwardRate2);
                }
                
                portfolioNpv+=temp.NPV();
            }
            
            
            if(ctptyDefault) {
                if(portfolioNpv>0) {
                    cva=(1-ctptyRecovery)*portfolioNpv/bankAccount;
                } else {
                    cva=0;
                }
            }
            
            if(invstDefault) {
                if(portfolioNpv<0) {
                    dva=-(1-invstRecovery)*portfolioNpv/bankAccount;
                } else {
                    dva=0;
                }
            }
            
            
        }
        
        if (!ctptyDefault) {
            cva=0;
        }
        if (!invstDefault) {
            dva=0;
        }
        
        cvaAccumulator.add(cva,pathWeight);
        dvaAccumulator.add(dva,pathWeight);
        
        cvaDvaSample[l][0]=cva;
        cvaDvaSample[l][1]=dva;
        cvaDvaSample[l][2]=x[0];
        cvaDvaSample[l][3]=x[1];
        
        
    }
    
    writeCsv(cvaDvaSample,"cvaDvaSampleSobol.csv");
    
    //reset evaluationdate
    Settings::instance().evaluationDate()=anchorDate;
    
    Matrix results(2,2);
    results[0][0]=cvaAccumulator.mean();
    results[1][0]=cvaAccumulator.errorEstimate();
    results[0][1]=dvaAccumulator.mean();
    results[1][1]=dvaAccumulator.errorEstimate();
    
    return results;
}




Matrix mcCvaSwapPortfolio(std::vector<ext::shared_ptr<VanillaSwap>> swapPortfolio,
                          ext::shared_ptr<MersenneTwisterUniformRng> uniformGenerator, Size nTimeSteps, Size nSamples,
                          Matrix rho, ext::shared_ptr<HullWhite> interestRate, ext::shared_ptr<CoxIngersollRoss> ctptyIntensity, ext::shared_ptr<CoxIngersollRoss> invstIntensity,
                          Real ctptyRecovery, Real invstRecovery) {
    
    Date anchorDate = Settings::instance().evaluationDate();
    
    Date longestMaturity =swapPortfolio[0]->maturityDate();
    for(int i=0;i<swapPortfolio.size();i++) {
        if(longestMaturity<swapPortfolio[i]->maturityDate()) {
            longestMaturity=swapPortfolio[i]->maturityDate();
        }
    }
    
    Actual365Fixed dummyDayCounter;
    TARGET dummyCalendar;
    
    TimeGrid timeGrid(dummyDayCounter.yearFraction(anchorDate,longestMaturity),nTimeSteps);
    
    InverseCumulativeRng<MersenneTwisterUniformRng,InverseCumulativeNormal> rng(*uniformGenerator);
    
    Matrix rhoChol = CholeskyDecomposition(rho);
    
    Natural nNormals = 3;//we need to have 4 normals for every step
    
    Real pathWeight=1;
    
    Statistics cvaAccumulator;
    Statistics dvaAccumulator;
    
    ext::shared_ptr<StochasticProcess> interestRateProcess = interestRate->dynamics()->process();
    ext::shared_ptr<StochasticProcess> ctptyIntensityProcess = ctptyIntensity->dynamics()->process();
    ext::shared_ptr<StochasticProcess> invstIntensityProcess = invstIntensity->dynamics()->process();
    
    Matrix cvaDvaSample(nSamples,3); //CVA, DVA and x
    
    for(int l=0;l<nSamples;l++) {
        Real ctptyU =1-uniformGenerator->nextReal();
        Real invstU =1-uniformGenerator->nextReal();
        
        Array x=interestRateProcess->initialValues();
        Array y=ctptyIntensityProcess->initialValues();
        Array z=invstIntensityProcess->initialValues();
        
        
        
        Real bankAccount=1;
        Real ctptySurvProb=1;
        Real invstSurvProb=1;
        
        bool ctptyDefault=false;
        bool invstDefault=false;
        
        Real portfolioNpv=0;
        Time tau;
        
        Real cva,dva;
        
        for(int i=0;i<nTimeSteps-1;i++) {//default on the last timestep does not count as a default
            Time t = timeGrid[i];
            Time dt =timeGrid.dt(i);
            Array dw(nNormals);
            for(int j=0;j<nNormals;j++) {
                dw[j]=rng.next().value;
            }
            Array dwCor(nNormals);
            dwCor=rhoChol*dw;
            
            Array tempX={dwCor[0]};
            Array tempY={dwCor[1]};
            Array tempZ={dwCor[2]};
            
            x=interestRateProcess->evolve(t,x,dt,tempX);
            y=ctptyIntensityProcess->evolve(t,y,dt,tempY);
            z=invstIntensityProcess->evolve(t,z,dt,tempZ);
            
            bankAccount/=interestRate->discountBond(t,t+dt,x);
            ctptySurvProb*=ctptyIntensity->discountBond(t,t+dt,y);
            invstSurvProb*=invstIntensity->discountBond(t,t+dt,z);
            
            if(ctptyU>=ctptySurvProb) {
                ctptyDefault=true;
                tau=t;
                break;
            }
            
            if(invstU>=invstSurvProb){
                invstDefault=true;
                tau=t;
                break;
            }
            
        }
        if(ctptyDefault||invstDefault) {
            Date defaultDate=anchorDate+tau*365;
            Settings::instance().evaluationDate()=defaultDate;
            
            std::vector<Date> dates;
            std::vector<Real> discFactors;
            
            for(int j=timeGrid.index(tau);j<nTimeSteps;j++) {
                dates.push_back(anchorDate+timeGrid[j]*365);
                discFactors.push_back(interestRate->discountBond(tau,timeGrid[j],x));
            }
            
            ext::shared_ptr<InterpolatedDiscountCurve<Cubic>> impliedTermStructure(
                                                                                   new InterpolatedDiscountCurve<Cubic>(dates,discFactors,dummyDayCounter,dummyCalendar));
            
            impliedTermStructure->enableExtrapolation();
            
            Handle<YieldTermStructure> impliedHandle(impliedTermStructure);
            for(int k=0;k<swapPortfolio.size();k++) {
                
                ext::shared_ptr<IborIndex> newIndex(
                                                    new IborIndex(swapPortfolio[k]->iborIndex()->familyName(),
                                                                  swapPortfolio[k]->iborIndex()->tenor(),
                                                                  swapPortfolio[k]->iborIndex()->fixingDays(),
                                                                  swapPortfolio[k]->iborIndex()->currency(),
                                                                  swapPortfolio[k]->iborIndex()->fixingCalendar(),
                                                                  swapPortfolio[k]->iborIndex()->businessDayConvention(),
                                                                  swapPortfolio[k]->iborIndex()->endOfMonth(),
                                                                  swapPortfolio[k]->iborIndex()->dayCounter(),
                                                                  impliedHandle));
                newIndex->clearFixings();//to remove fixings from other paths
                VanillaSwap temp(swapPortfolio[k]->type(),
                                 swapPortfolio[k]->nominal(),
                                 swapPortfolio[k]->fixedSchedule(),
                                 swapPortfolio[k]->fixedRate(),
                                 swapPortfolio[k]->fixedDayCount(),
                                 swapPortfolio[k]->floatingSchedule(),
                                 newIndex,
                                 swapPortfolio[k]->spread(),
                                 swapPortfolio[k]->floatingDayCount(),
                                 swapPortfolio[k]->paymentConvention());
                //discounting using the forward curve
                Handle<YieldTermStructure> disc =newIndex->forwardingTermStructure();
                ext::shared_ptr<PricingEngine> engine(new
                                                      DiscountingSwapEngine(disc));
                
                temp.setPricingEngine(engine);
                
                Date dateBefore =temp.floatingSchedule().previousDate(defaultDate);
                Date dateAfter1 =temp.floatingSchedule().nextDate(defaultDate);
                Date dateAfter2 =temp.floatingSchedule().nextDate(dateAfter1);
                
                //approximate the last fixing as the current simple forward rate
                Rate forwardRate1 = impliedTermStructure->forwardRate(defaultDate,dateAfter1, temp.floatingDayCount(), Simple);
                Rate forwardRate2 = impliedTermStructure->forwardRate(dateAfter1,dateAfter2, temp.floatingDayCount(), Simple);
                
                if(dateBefore!=Date()) {//if not null date
                    Date fixing1 =newIndex->fixingDate(dateBefore);
                    newIndex->addFixing(fixing1,forwardRate1);
                }
                
                Date fixing2 =newIndex->fixingDate(dateAfter1);
                if(fixing2<defaultDate) {
                    newIndex->addFixing(fixing2,forwardRate2);
                }
                
                portfolioNpv+=temp.NPV();
            }
            
            
            if(ctptyDefault) {
                if(portfolioNpv>0) {
                    cva=(1-ctptyRecovery)*portfolioNpv/bankAccount;
                } else {
                    cva=0;
                }
            }
            
            if(invstDefault) {
                if(portfolioNpv<0) {
                    dva=-(1-invstRecovery)*portfolioNpv/bankAccount;
                } else {
                    dva=0;
                }
            }
            
            
        }
        
        if (!ctptyDefault) {
            cva=0;
        }
        if (!invstDefault) {
            dva=0;
        }
        
        cvaAccumulator.add(cva,pathWeight);
        dvaAccumulator.add(dva,pathWeight);
        
        cvaDvaSample[l][0]=cva;
        cvaDvaSample[l][1]=dva;
        cvaDvaSample[l][2]=x[0];
        
        
    }
    
    writeCsv(cvaDvaSample,"cvaDvaSampleHw.csv");
    
    //reset evaluationdate
    Settings::instance().evaluationDate()=anchorDate;
    
    Matrix results(2,2);
    results[0][0]=cvaAccumulator.mean();
    results[1][0]=cvaAccumulator.errorEstimate();
    results[0][1]=dvaAccumulator.mean();
    results[1][1]=dvaAccumulator.errorEstimate();
    
    return results;
}


