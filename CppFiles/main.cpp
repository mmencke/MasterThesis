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

int main() {
    //std::locale::global(std::locale("en_US.UTF-8")); //sets global number formatting to US
    
    bool run_section = false;
    if(run_section) {
        newSection("Check Version");
        checkVersion();
    }
    
    /*******************************************************************/
    /* CORRELATED NORMALS                                              */
    /*******************************************************************/
    
    run_section = true;
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
    /* DISCRETISATION OF CIR                                       */
    /*******************************************************************/
    
    run_section = true;
    if(run_section) {
        newSection("Discretisation of CIR");
        discretiseCir();
    }
    
    
    /*******************************************************************/
    /* CALIBRATE CIR TO OLD MARKET DATA*/
    /*******************************************************************/
    
    run_section = false;
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
    /* MONTE CARLO CVA AND DVA                          */
    /*******************************************************************/
    
    run_section=true;
    if(run_section) {
        newSection("CVA and DVA by Monte Carlo");
        mcCva();
    }
    
    
    
    return 0;
}
