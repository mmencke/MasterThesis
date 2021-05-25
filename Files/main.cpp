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
#include "DiscretiseVasicek.hpp"
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
    /* DISCRETISATION OF VASICEK                                       */
    /*******************************************************************/
    
    run_section = false;
    if(run_section) {
        newSection("Discretisation of Vasicek");
        discretiseVasicek();
    }
    
    /*******************************************************************/
    /* CALIBRATE CIR TO OLD MARKET DATA*/
    /*******************************************************************/
    
    run_section = true;
    if(run_section) {
        newSection("Calibrate CIR to Old Market Data");
        calibrateCir();
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
