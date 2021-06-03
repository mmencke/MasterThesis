//
//  CorrelatedNormals.hpp
//  Master Thesis
//
//  Created by Magnus Mencke on 26/05/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef CorrelatedNormals_hpp
#define CorrelatedNormals_hpp

#include "headers.hpp"

inline void correlatedNormals() {
    

    unsigned long seed =1;

    MersenneTwisterUniformRng uniformGenerator(seed);

    InverseCumulativeRng<MersenneTwisterUniformRng,InverseCumulativeNormal> rng(uniformGenerator);

    Natural n = 10000;

    Matrix rho(3,3);
    rho[0][0]=rho[1][1]=rho[2][2]=1;
    rho[0][1]=rho[1][0]=0;
    rho[0][2]=rho[2][0]=0.5;
    rho[1][2]=rho[2][1]=-0.5;

    std::cout << std::endl;
    printMatrix(rho);
    std::cout << std::endl;

    Matrix rhoChol=CholeskyDecomposition(rho);

    std::cout << std::endl;
    printMatrix(rhoChol);
    std::cout << std::endl;

    Matrix correlatedNormals(n,3);

    for(int k=0;k<n;k++) {
        Array dw(3);
        for(int j=0;j<3;j++) {
            dw[j]=rng.next().value;
        }
        Array dwCor(3);
        dwCor=rhoChol*dw;
        
        for(int j=0;j<3;j++) {
            correlatedNormals[k][j]=dwCor[j];
        }
    }

    writeCsv(correlatedNormals,"correlatedNormals.csv");
}
#endif /* CorrelatedNormals_hpp */
