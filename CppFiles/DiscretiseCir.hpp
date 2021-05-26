//
//  DiscretiseCir.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 25/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef DiscretiseCir_hpp
#define DiscretiseCir_hpp

#include "headers.hpp"
#include "ExtraCirProcesses.hpp"

inline void discretiseCir() {
    BigNatural seed = 1;
    Size nTimeSteps = 252;
    Size nSamples = 10000;
    Statistics statistics;
    
    bool brownianBridge = false;
    bool antitheticVariate = false;
    bool withFellerConstraint = true;
    
    //Parameters from [Brigo, Morini and Pallavicini (2013), p. 125]
    Real kappa=0.4;
    Real mu = 0.026;
    Real v = 0.14;
    Real y0 = 0.0165;
    
    TimeGrid timeGrid(1,nTimeSteps);
    
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(nTimeSteps, seed);
    
    Matrix cirDistribution(nSamples,5);
    
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    
    /*******************************************************************/
    /* SIMPLE EULER DISCRETISATION                                     */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    
    ext::shared_ptr<StochasticProcess1D> cirProcessSimple(new CirProcess(kappa, v, y0, mu, CirProcess::None));
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cirRsgSimple(new SingleVariate<PseudoRandom>::path_generator_type(cirProcessSimple,timeGrid, rsg, brownianBridge));
    
    
    for(int i=0;i<nSamples;i++) {
        Path path = (cirRsgSimple->next()).value;
        cirDistribution[i][0] =path[nTimeSteps];
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << boost::format("%-25s %-12s %-8s %-3s") % "Simple Euler:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
                                                          
    /*******************************************************************/
    /* FULL TRUNCATION SCHEME                                          */
    /*******************************************************************/
                  
    start = std::chrono::steady_clock::now();
    
    ext::shared_ptr<StochasticProcess1D> cirProcessTrunc(new CirProcess(kappa, v, y0, mu, CirProcess::FullTruncation));
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cirRsgTrunc(new SingleVariate<PseudoRandom>::path_generator_type(cirProcessTrunc,timeGrid, rsg, brownianBridge));
    
    
    for(int i=0;i<nSamples;i++) {
        Path path = (cirRsgTrunc->next()).value;
        cirDistribution[i][1] =path[nTimeSteps];
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << boost::format("%-25s %-12s %-8s %-3s") % "Full Truncation:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
    /*******************************************************************/
    /* SIMULATE SQUARE ROOT                                            */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    
    ext::shared_ptr<StochasticProcess1D> cirProcessSqrt(new CirHelperProcess(kappa, v, y0, mu));
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cirRsgSqrt(new SingleVariate<PseudoRandom>::path_generator_type(cirProcessSqrt,timeGrid, rsg, brownianBridge));

    for(int i=0;i<nSamples;i++) {
        Path path = (cirRsgSqrt->next()).value;
        cirDistribution[i][2] =std::pow(path[nTimeSteps],2);
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << boost::format("%-25s %-12s %-8s %-3s") % "Simulate Square Root:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
    /*******************************************************************/
    /* QUADRATIC EXPONENTIAL SCHEME                                    */
    /*******************************************************************/
    
    start = std::chrono::steady_clock::now();
    
    ext::shared_ptr<StochasticProcess1D> cirProcessQe(new CoxIngersollRossProcess(kappa, v, y0, mu));
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cirRsgQe(new SingleVariate<PseudoRandom>::path_generator_type(cirProcessQe,timeGrid, rsg, brownianBridge));
    
    for(int i=0;i<nSamples;i++) {
        Path path = (cirRsgQe->next()).value;
        cirDistribution[i][3] =path[nTimeSteps];
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << boost::format("%-25s %-12s %-8s %-3s") % "Quadratic Exponential:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
    /*******************************************************************/
    /* EXACT SIMULATION                                                */
    /*******************************************************************/
    
    
    start = std::chrono::steady_clock::now();
    
    ext::shared_ptr<StochasticProcess1D> cirProcessExact(new CirProcess(kappa, v, y0, mu, CirProcess::Exact));
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> cirRsgExact(new SingleVariate<PseudoRandom>::path_generator_type(cirProcessExact,timeGrid, rsg, brownianBridge));
    
    
    for(int i=0;i<nSamples;i++) {
        Path path = (cirRsgExact->next()).value;
        cirDistribution[i][4] =path[nTimeSteps];
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout << boost::format("%-25s %-12s %-8s %-3s") % "Exact Simulation:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
    
    writeCsv(cirDistribution,"cirDistribution.csv");
    
    /*******************************************************************/
    /* QUADRATIC EXPONENTIAL TEST OF CONVERGENCE                       */
    /*******************************************************************/
    
    BigNatural nSamplesQe = 1000000;
    
    Matrix cirQeConvergence(nSamplesQe,1);
    
    start = std::chrono::steady_clock::now();
    
    for(int i=0;i<nSamplesQe;i++) {
        Path path = (cirRsgQe->next()).value;
        cirQeConvergence[i][0] =path[nTimeSteps];
    }
    
    end = std::chrono::steady_clock::now();
    
    diff = end - start;
    
    std::cout  << std::endl << boost::format("%-25s %-12s %-8s %-3s") %  "QE Convergence:" % "Run time =" % std::chrono::duration <double, std::milli> (diff).count() % "ms" << std::endl;
    
    
    writeCsv(cirQeConvergence,"cirQeConvergence.csv");
    
}


#endif /* DiscretiseCir_hpp */
