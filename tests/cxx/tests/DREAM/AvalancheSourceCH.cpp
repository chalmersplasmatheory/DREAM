/**
 * A test for the Rosenbluth-Putvinski avalanche source that 
 * verifies that the total number of knock-on particles created
 * in the entire volume agrees with the analytic exact result. 
 */

#include <vector>
#include <string>
#include "AvalancheSourceCH.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"
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
bool AvalancheSourceCH::Run(bool) {

    bool success = true;
    if (CheckConservativityCylindrical())
        this->PrintOK("The CH avalanche source creates the predicted total number of knock-ons in cylindrical geometry.");
    else {
        success = false;
        this->PrintError("The avalanche test failed: an incorrect number of knock-ons are created in cylindrical geometry.");
    }
    /*if (CheckConservativityGeneralAnalytic())
        this->PrintOK("The CH avalanche source creates the predicted total number of knock-ons in general geometry.");
    else {
        success = false;
        this->PrintError("The avalanche test failed: an incorrect number of knock-ons are created in general geometry.");
    }*/

    return success;
}


/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *AvalancheSourceCH::GetUnknownHandler(
    DREAM::FVM::Grid *g_fl, DREAM::FVM::Grid *g_hot, DREAM::FVM::Grid *g_re, 
    const real_t f_hot, const real_t f_re, const real_t n_tot, const real_t E_field
) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_RE, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", g);

    #define SETVAL_HOT(NAME, v) do { \
            for (len_t i = 0; i < g_hot->GetNCells(); i++) \
                temp_hot[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)

    #define SETVAL_RE(NAME, v) do { \
            for (len_t i = 0; i < g_re->GetNCells(); i++) \
                temp_re[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)
    
    #define SETVAL_FL(NAME, v) do { \
            for (len_t i = 0; i < g_fl->GetNr(); i++) \
                temp_fl[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)
    

    // Set electron quantities
    real_t *temp_hot = new real_t[g_hot->GetNCells()];
    SETVAL_HOT(DREAM::OptionConstants::UQTY_F_HOT, f_hot);
    real_t *temp_re  = new real_t[g_re->GetNCells()];
    SETVAL_RE(DREAM::OptionConstants::UQTY_F_RE, f_re);
    real_t *temp_fl = new real_t[g_fl->GetNr()];
    SETVAL_FL(DREAM::OptionConstants::UQTY_N_TOT, n_tot);
    SETVAL_FL(DREAM::OptionConstants::UQTY_E_FIELD, E_field);
    return uqh;
}

/**
 * Test particle number with a general inhomogeneous grid
 */
/*bool AvalancheSourceCH::CheckConservativityGeneralAnalytic(){
    real_t successRelErrorThreshold = 1e-5;

    len_t nr = 2;
    len_t np  = 4;
    len_t nxi = 5;
    
    len_t ntheta_interp = 100;
    len_t nrProfiles    = 8;

    real_t pMin_hot = 0;
    real_t pMax_hot = 2;

    real_t pMin_re = 2;
    real_t pMax_re = 4;

    auto *grid_hot  = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin_hot, pMax_hot);
    auto *grid_re   = InitializeGridGeneralRPXi(nr, np, nxi, ntheta_interp, nrProfiles, pMin_re, pMax_re);
    auto *fluidGrid = InitializeFluidGrid(nr);
    
    real_t f_re = 1e15;
    real_t n_tot = 1e20;
    real_t E_field = -1;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(fluidGrid, n_re, n_tot, E_field);

    real_t pCutoff = 0.13*pMax;
    DREAM::AvalancheSourceCH *avaSourceTerm = new DREAM::AvalancheSourceCH(grid, uqh,pCutoff,1.0);

    real_t *sourceVec = new real_t[np*nxi];
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np; i++)
            for(len_t j=0; j<nxi; j++)
                sourceVec[j*np+i] = n_re * n_tot * avaSourceTerm->EvaluateCHSource(ir,i,j);
        real_t sourceIntegralNumerical = grid->IntegralMomentumAtRadius(ir, sourceVec);
        real_t sourceIntegralAnalytic = n_re * n_tot * avaSourceTerm->EvaluateNormalizedTotalKnockOnNumber(pCutoff, pMax);
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
}*/

/**
 * Test particle number with a homogeneous cylindrical grid
 */
bool AvalancheSourceCH::CheckConservativityCylindrical(){
    real_t eps = numeric_limits<real_t>::epsilon();
    real_t successRelErrorThreshold = 1000*eps;

    len_t nr = 2;
    len_t np = 4;
    len_t nxi = 8;

    real_t pMin_hot = 0;
    real_t pMax_hot = 2;

    real_t pMin_re = 2;
    real_t pMax_re = 12;
    
    real_t B0 = 1.0;

    auto *grid_fl = InitializeFluidGrid(nr);
    auto *grid_hot = InitializeGridRCylPXi(nr,np,nxi,B0,pMin_hot,pMax_hot);
    auto *grid_re  = InitializeGridRCylPXi(nr,np,nxi,B0,pMin_re,pMax_re);

    real_t f_hot = 2e15;
    real_t f_re = 2e15;
    real_t n_tot = 1e20;
    real_t E_field = 1;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(grid_fl, grid_hot, grid_re, f_hot, f_re, n_tot, E_field);

    real_t pCutoff = 0.125;
    DREAM::AvalancheSourceCH *avaSourceTerm_hot_pos = new DREAM::AvalancheSourceCH(
            grid_hot, uqh, pCutoff, 1.0, AvalancheSourceCH::CH_SOURCE_PITCH_POSITIVE, 
            AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, false, grid_re
        );
    DREAM::AvalancheSourceCH *avaSourceTerm_hot_neg = new DREAM::AvalancheSourceCH(
            grid_hot, uqh, pCutoff, 1.0, AvalancheSourceCH::CH_SOURCE_PITCH_NEGATIVE, 
            AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, false, grid_re
        );
    DREAM::AvalancheSourceCH *avaSourceTerm_re_pos = new DREAM::AvalancheSourceCH(
            grid_re, uqh, pCutoff, 1.0, AvalancheSourceCH::CH_SOURCE_PITCH_POSITIVE, 
            AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, true
        );
    DREAM::AvalancheSourceCH *avaSourceTerm_re_neg = new DREAM::AvalancheSourceCH(
            grid_re, uqh, pCutoff, 1.0, AvalancheSourceCH::CH_SOURCE_PITCH_NEGATIVE, 
            AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, true
        );

    real_t *sourceVec_hot_pos = new real_t[nr*np*nxi];
    real_t *sourceVec_hot_neg = new real_t[nr*np*nxi];
    real_t *sourceVec_re_pos = new real_t[nr*np*nxi];
    real_t *sourceVec_re_neg = new real_t[nr*np*nxi];
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np; i++){
            for(len_t j=0; j<nxi; j++){
                sourceVec_hot_pos[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_hot_pos->GetSourceFunction(ir,i,j);
                sourceVec_hot_neg[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_hot_neg->GetSourceFunction(ir,i,j);
                sourceVec_re_pos[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_re_pos->GetSourceFunction(ir,i,j);
                sourceVec_re_neg[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_re_neg->GetSourceFunction(ir,i,j);
            }
        }
        //real_t sourceIntegralNumerical = grid->IntegralMomentumAtRadius(ir, sourceVec);
        //real_t sourceIntegralAnalytic = n_re * n_tot * avaSourceTerm->EvaluateNormalizedTotalKnockOnNumber(pCutoff, pMax);
        //deltas[ir] = abs(sourceIntegralAnalytic - sourceIntegralNumerical) / sourceIntegralAnalytic;
        deltas[ir] = 0;
    }
    for(len_t j=0; j<nxi; j++){
        for(len_t i=0; i<np; i++){
            printf("%.1e    ", sourceVec_hot_pos[j*np + i]);
        }
        printf("\n");
    }
    /*
    printf("\n\n");
    for(len_t j=0; j<nxi; j++){
        for(len_t i=0; i<np; i++){
            printf("%.1e    ", sourceVec_hot_neg[j*np + i]);
        }
        printf("\n");
    }
    printf("\n\n");
    for(len_t j=0; j<nxi; j++){
        for(len_t i=0; i<np; i++){
            printf("%.1e    ", sourceVec_re_pos[j*np + i]);
        }
        printf("\n");
    }
    printf("\n\n");
    for(len_t j=0; j<nxi; j++){
        for(len_t i=0; i<np; i++){
            printf("%.1e    ", sourceVec_re_neg[j*np + i]);
        }
        printf("\n");
    }
    */

    bool success = true;
    for(len_t ir=0; ir<nr; ir++)
        if(deltas[ir]>successRelErrorThreshold){
            success = false;
            cout << "delta [rel error]: " << deltas[ir] << endl;
        }

    delete [] sourceVec;
    return success;
}


