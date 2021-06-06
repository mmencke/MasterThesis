//
//  main.cpp
//  Master Thesis
//
//  Created by Magnus Mencke on 29/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#include "utility.hpp"

#include "CheckVersion.hpp"
#include "CorrelatedNormals.hpp"
#include "DiscretiseVasicek.hpp"
#include "DiscretiseCir.hpp"
#include "CalibrateCir.hpp"
#include "Calibration.hpp"
#include "MonteCarlo.hpp"
#include "SwapNpvShift.hpp"
#include "MonteCarloSample.hpp"

int main() {
    
    bool run_section;
    
    run_section = false;
    if(run_section) {
        newSection("Check Version");
        checkVersion();
    }
    
    /*******************************************************************/
    /* CORRELATED NORMALS                                              */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Correlated Normals");
        correlatedNormals();
    }
    
    /*******************************************************************/
    /* DISCRETISATION OF VASICEK                                       */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Discretisation of Vasicek");
        discretiseVasicek();
    }
    
    /*******************************************************************/
    /* DISCRETISATION OF CIR                                           */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Discretisation of CIR");
        discretiseCir();
    }
    
    /*******************************************************************/
    /* CALIBRATE CIR TO OLD MARKET DATA                                */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Calibrate CIR to Old Market Data");
        calibrateCir();
    }
    
    /*******************************************************************/
    /* CALIBRATION TO NEW MARKET DATA                                  */
    /*******************************************************************/

    run_section = false;
    if(run_section) {
        newSection("Calibration of Interest Rate and Intensity");
        calibration();
    }
    
    /*******************************************************************/
    /* MONTE CARLO CVA AND DVA                                         */
    /*******************************************************************/
    
    run_section = true;
    if(run_section) {
        newSection("CVA and DVA by Monte Carlo");
        mcCva();
    }
    
    /*******************************************************************/
    /* SWAP NPV SHIFT                                                  */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Changes in Swap NPV by shifts in the Gaussian Factors");
        swapNpvShift();
    }
    
    /*******************************************************************/
    /* MONTE CARLO SAMPLE                                              */
    /*******************************************************************/
    
    // Every time the Monte Carlo simulation is run, the samples are output along with the realisation of the Gaussian Factors.
    // We want to investigate the 30y independent ATM case, so this is run here alone to overwrite the csv.
    
    run_section = false;
    if(run_section) {
        newSection("Monte Carlo Sample (Single Maturity)");
        monteCarloSample();
    }
    return 0;
}
