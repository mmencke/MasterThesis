//
//  SimulateVasicek.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef SimulateVasicek_hpp
#define SimulateVasicek_hpp



#include <iostream>
#include <ql/models/shortrate/onefactormodels/vasicek.hpp>
#include <ql/methods/montecarlo/montecarlomodel.hpp>

#include "utility.hpp"

using namespace QuantLib;

inline void simulateVasicek() {
    Real kappa = 0.05;
    Real theta =0.03;
    Real sigma = 0.01; //absolute volatility
    Real r0 = 0.03;
    
    BigNatural seed = 1;
    Size n_time_steps = 250;
    Size n_samples = 100000;
    Statistics statistics_accumulator;
    
    TimeGrid time_grid(1,n_time_steps);
    
    bool brownian_bridge = false;
    bool antithetic_variate = false;
    
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(n_time_steps, seed);
    
    Matrix vasicek_distribution(n_samples,1);
    
    auto start = std::chrono::steady_clock::now();
    
    Vasicek vasicek_model(r0, kappa, theta,sigma);
    
    ext::shared_ptr<StochasticProcess1D> vasicek_process_ptr = vasicek_model.dynamics()->process();
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> vasicek_rsg_ptr(new SingleVariate<PseudoRandom>::path_generator_type(vasicek_process_ptr,time_grid, rsg, brownian_bridge));
    
    for(int i=0;i<n_samples;i++) {
        Path path = (vasicek_rsg_ptr->next()).value;
        vasicek_distribution[i][0] =path[n_time_steps-1];
    }
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "Simple Euler discretization: run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    
    writeCsv(vasicek_distribution,"vasicek_dist.csv");
}

#endif /* SimulateVasicek_hpp */
