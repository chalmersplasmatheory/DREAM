/**
 * A test for the Rosenbluth-Putvinski avalanche source that 
 * verifies that the total number of knock-on particles created
 * in the entire volume agrees with the analytic exact result. 
 */

#include <vector>
#include <string>
#include "AvalancheSourceRP.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>

using namespace DREAMTESTS::_DREAM;
using namespace std;




/**
 * Run this test.
 */
bool AvalancheSourceRP::Run(bool) {

    bool success = true;
    if (CheckConservativity())
        this->PrintOK("The RP avalanche source creates the exact analytic total number of knock-ons.");
    else {
        success = false;
        this->PrintError("The avalanche test failed: an incorrect number of knock-ons are created.");
    }

    return success;
}


/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *AvalancheSourceRP::GetUnknownHandler(DREAM::FVM::Grid *g, 
    const real_t n_re, const real_t n_tot) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_RE, g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, g);

    #define SETVAL(NAME, v) do { \
            for (len_t i = 0; i < g->GetNr(); i++) \
                temp[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)

    // Set electron quantities
    real_t *temp = new real_t[g->GetNr()];
    SETVAL(DREAM::OptionConstants::UQTY_N_RE, n_re);
    SETVAL(DREAM::OptionConstants::UQTY_N_TOT, n_tot);
    
    return uqh;
}

bool AvalancheSourceRP::CheckConservativity(){
    real_t successRelErrorThreshold = 1e-3;

    len_t nr = 2;
    len_t np = 3;
    len_t nxi = 6;

    len_t ntheta_ref = 12;
    len_t ntheta_interp = 8;
    len_t nrProfiles = 5;

    real_t pMin = 0;
    real_t pMax = 0.4;

    auto *grid = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_ref, ntheta_interp, nrProfiles, pMin, pMax);
    auto *fluidGrid = InitializeFluidGrid(nr);
    real_t n_re = 1e13;
    real_t n_tot = 1e20;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(fluidGrid, n_re, n_tot);

    real_t pCutoff = 0.13*pMax;
    DREAM::AvalancheSourceRP *avaSourceTerm = new DREAM::AvalancheSourceRP(grid, uqh,pCutoff);

    real_t **sourceVec = new real_t*[nr];
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        sourceVec[ir] = new real_t[np*nxi];
        for(len_t i=0; i<np; i++){
            for(len_t j=0; j<nxi; j++){
                sourceVec[ir][j*np+i] = n_re * n_tot * avaSourceTerm->EvaluateRPSource(ir,i,j);
                if( sourceVec[ir][j*np+i] ){
                    //cout << "Non-zero contribution at ir = " << ir << ", i = " << i << ", j = " << j << "." << endl;
                    break;
                }
            }
        }
        real_t sourceIntegralNumerical = grid->IntegralMomentumAtRadius(ir, sourceVec[ir]);
        real_t sourceIntegralAnalytic = avaSourceTerm->EvaluateTotalKnockOnNumber(ir, pCutoff, pMax);
        deltas[ir] = abs(sourceIntegralAnalytic - sourceIntegralNumerical) / sourceIntegralAnalytic;
/*
        cout << "Analytic: " << sourceIntegralAnalytic << endl;
        cout << "Numeric: " << sourceIntegralNumerical << endl;
        cout << "Ratio: " << sourceIntegralAnalytic / sourceIntegralNumerical << endl;
*/
        delete [] sourceVec[ir];
    }

    bool success = true;
    for(len_t ir=0; ir<nr; ir++)
        if(deltas[ir]>successRelErrorThreshold){
            success = false;
            cout << "delta [rel error]: " << deltas[ir] << endl;
        }

    delete [] sourceVec;
    return success;

}


