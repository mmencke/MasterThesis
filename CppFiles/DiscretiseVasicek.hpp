//
//  SimulateVasicek.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef DiscretiseVasicek_hpp
#define DiscretiseVasicek_hpp


#include "headers.hpp"



inline void discretiseVasicek() {
    Real kappa = 0.05;
    Real theta =0.03;
    Real sigma = 0.01; //absolute volatility
    Real r0 = 0.03;
    
    BigNatural seed = 1;
    Size nTimeSteps = 250;
    Size nSamples = 100000;
    Statistics statistics;
    
    TimeGrid timeGrid(1, nTimeSteps);
    
    bool brownianBridge = false;
    bool antitheticVariate = false;
    
    PseudoRandom::rsg_type rsg = PseudoRandom::make_sequence_generator(nTimeSteps, seed);
    
    Matrix vasicekDistribution(nSamples,1);
    
    auto start = std::chrono::steady_clock::now();
    
    Vasicek vasicekModel(r0, kappa, theta, sigma);
    
    ext::shared_ptr<StochasticProcess1D> vasicekProcess = vasicekModel.dynamics()->process();
    
    ext::shared_ptr<SingleVariate<PseudoRandom>::path_generator_type> vasicekRsg(new SingleVariate<PseudoRandom>::path_generator_type(vasicekProcess, timeGrid, rsg, brownianBridge));
    
    for(int i=0;i<nSamples;i++) {
        Path path = (vasicekRsg->next()).value;
        vasicekDistribution[i][0] =path[nTimeSteps];
    }
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "Run time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    
    writeCsv(vasicekDistribution,"vasicekDistribution.csv");
}

#endif /* SimulateVasicek_hpp */
