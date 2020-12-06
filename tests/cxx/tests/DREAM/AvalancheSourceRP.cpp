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
    if (CheckConservativityCylindrical())
        this->PrintOK("The RP avalanche source creates the predicted total number of knock-ons in cylindrical geometry.");
    else {
        success = false;
        this->PrintError("The avalanche test failed: an incorrect number of knock-ons are created in cylindrical geometry.");
    }
    if (CheckConservativityGeneralAnalytic())
        this->PrintOK("The RP avalanche source creates the predicted total number of knock-ons in general geometry.");
    else {
        success = false;
        this->PrintError("The avalanche test failed: an incorrect number of knock-ons are created in general geometry.");
    }

    return success;
}


/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *AvalancheSourceRP::GetUnknownHandler(DREAM::FVM::Grid *g, 
    const real_t n_re, const real_t n_tot, const real_t E_field) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_RE, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", g);

    #define SETVAL(NAME, v) do { \
            for (len_t i = 0; i < g->GetNr(); i++) \
                temp[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)

    // Set electron quantities
    real_t *temp = new real_t[g->GetNr()];
    SETVAL(DREAM::OptionConstants::UQTY_N_RE, n_re);
    SETVAL(DREAM::OptionConstants::UQTY_N_TOT, n_tot);
    SETVAL(DREAM::OptionConstants::UQTY_E_FIELD, E_field);
    return uqh;
}

/**
 * Test particle number with a general inhomogeneous grid
 */
bool AvalancheSourceRP::CheckConservativityGeneralAnalytic(){
    real_t successRelErrorThreshold = 1e-5;

    len_t nr = 2;
    len_t np = 4;
    len_t nxi = 5;

    len_t ntheta_interp = 100;
    len_t nrProfiles = 8;

    real_t pMin = 0;
    real_t pMax = 2;

    auto *grid = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin, pMax);
    auto *fluidGrid = InitializeFluidGrid(nr);
    real_t n_re = 1e13;
    real_t n_tot = 1e20;
    real_t E_field = -1;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(fluidGrid, n_re, n_tot, E_field);

    real_t pCutoff = 0.13*pMax;
    DREAM::AvalancheSourceRP *avaSourceTerm = new DREAM::AvalancheSourceRP(grid, uqh,pCutoff,1.0);

    real_t *sourceVec = new real_t[np*nxi];
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np; i++)
            for(len_t j=0; j<nxi; j++)
                sourceVec[j*np+i] = n_re * n_tot * avaSourceTerm->EvaluateRPSource(ir,i,j);
        real_t sourceIntegralNumerical = grid->IntegralMomentumAtRadius(ir, sourceVec);
        real_t sourceIntegralAnalytic = n_re * n_tot * avaSourceTerm->EvaluateNormalizedTotalKnockOnNumber(grid->GetRadialGrid()->GetFSA_B(ir), pCutoff, pMax);
        deltas[ir] = abs(sourceIntegralAnalytic - sourceIntegralNumerical) / sourceIntegralAnalytic;
        //cout << "numerical: " << sourceIntegralNumerical << endl;
        //cout << "analytic: " << sourceIntegralAnalytic << endl;
    }

    bool success = true;
    for(len_t ir=0; ir<nr; ir++){
        if(deltas[ir]>successRelErrorThreshold){
            success = false;
            cout << "delta [rel error]: " << deltas[ir] << endl;
        }
    }

    delete [] sourceVec;
    return success;
}

/**
 * Test particle number with a homogeneous cylindrical grid
 */
bool AvalancheSourceRP::CheckConservativityCylindrical(){
    real_t eps = numeric_limits<real_t>::epsilon();
    real_t successRelErrorThreshold = 1000*eps;

    len_t nr = 2;
    len_t np = 3;
    len_t nxi = 6;

    real_t pMin = 0;
    real_t pMax = 0.4;
    
    real_t B0 = 1.0;

    auto *grid = InitializeGridRCylPXi(nr,np,nxi,B0,pMin,pMax);

    auto *fluidGrid = InitializeFluidGrid(nr);
    real_t n_re = 1e13;
    real_t n_tot = 1e20;
    real_t E_field = 1;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(fluidGrid, n_re, n_tot, E_field);

    real_t pCutoff = 0.13*pMax;
    DREAM::AvalancheSourceRP *avaSourceTerm = new DREAM::AvalancheSourceRP(grid, uqh,pCutoff,1.0);

    real_t *sourceVec = new real_t[np*nxi];
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np; i++)
            for(len_t j=0; j<nxi; j++)
                sourceVec[j*np+i] = n_re * n_tot * avaSourceTerm->EvaluateRPSource(ir,i,j);                
        real_t sourceIntegralNumerical = grid->IntegralMomentumAtRadius(ir, sourceVec);
        real_t sourceIntegralAnalytic = n_re * n_tot * avaSourceTerm->EvaluateNormalizedTotalKnockOnNumber(grid->GetRadialGrid()->GetFSA_B(ir), pCutoff, pMax);
        deltas[ir] = abs(sourceIntegralAnalytic - sourceIntegralNumerical) / sourceIntegralAnalytic;
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


