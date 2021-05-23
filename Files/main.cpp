//
//  main.cpp
//  Master Thesis
//
//  Created by Magnus Mencke on 29/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "checkversion.hpp"
#include "utility.hpp"
#include "SimpleSwapPricing.hpp"
#include "SimulateVasicek.hpp"
#include "VasicekSwapNoCr.hpp"
#include "VasicekSwapCr.hpp"
#include "CalibrateG2.hpp"
#include "cdsStuff.hpp"
#include "McCirCds.hpp"
#include "CalibrateCir.hpp"
#include "ApproxZcb.hpp"
#include "MonteCarlo.hpp"
#include "Calibration.hpp"

#include <random>



using namespace QuantLib;

int main() {
    //std::locale::global(std::locale("en_US.UTF-8")); //sets global number formatting to US
    
    bool run_section = false;
    if(run_section) {
        newSection("Check Version");
        checkVersion();
    }
    
    /*******************************************************************/
    /* SIMPLE SWAP PRICING                                             */
    /*******************************************************************/
    
    run_section=false;
    if(run_section) {
        newSection("Simple Swap Pricing");
        simpleSwapPricing();
    }
    
    /*******************************************************************/
    /* SIMULATE FROM VASICEK                                               */
    /*******************************************************************/
    run_section=false;
    if(run_section) {
        newSection("Simulate from Vasicek");
        simulateVasicek();
    }
    /*******************************************************************/
    /* VASICEK NO CR                                               */
    /*******************************************************************/
    run_section=false;
    if(run_section) {
        newSection("Swap in Vasicek with No Credit Risk");
        vasicekSwapNoCr();
    }
    
    /*******************************************************************/
    /* VASICEK CR                              */
    /*******************************************************************/
    
    run_section=false;
    if(run_section) {
        newSection("Swap in Vasicek with Credit Risk");
        vasicekSwapCr();
    }
    
    /*******************************************************************/
    /* CALIBRATE G2                              */
    /*******************************************************************/

    run_section=false;
    if(run_section) {
        newSection("Calibration of G2");
        calibrateG2();
    }
    
    
    run_section=false;
    if(run_section) {
        newSection("CDS Stuff");
        cdsStuff();
    }
    
    run_section=false;
    if(run_section) {
        newSection("CDS in CIR by MC");
        mcCirCds();
    }
    
    /*******************************************************************/
    /* Calibration of Interest Rate and Intensity                      */
    /*******************************************************************/

    run_section=false;
    if(run_section) {
        newSection("Calibration of Interest Rate and Intensity");
        calibration();
    }
    
    /*******************************************************************/
    /* CALIBRATE G2                              */
    /*******************************************************************/
    
    run_section=true;
    if(run_section) {
        newSection("Calibrate CIR");
        calibrateCir();
    }
    /*******************************************************************/
    /* Approximation of Bank Account / ZCB                           */
    /*******************************************************************/
    
    run_section=false;
    if(run_section) {
        newSection("Approx. ZCB");
        approxZcb();
    }
    
    /*******************************************************************/
    /* MONTE CARLO CVA                          */
    /*******************************************************************/
    
    run_section=false;
    if(run_section) {
        newSection("Monte Carlo Cva");
        mcCva();
    }
    
    
    
    return 0;
}
