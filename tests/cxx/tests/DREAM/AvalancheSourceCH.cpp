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
#include "DREAM/Constants.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Grid/Grid.hpp"

#include <iostream>

using namespace DREAMTESTS::_DREAM;
using namespace std;




/**
 * Run this test.
 */
bool AvalancheSourceCH::Run(bool) {
    gsl_ws1 = gsl_integration_workspace_alloc(1000);
    gsl_ws2 = gsl_integration_workspace_alloc(1000);
    gsl_ws3 = gsl_integration_workspace_alloc(1000);

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

    gsl_integration_workspace_free(gsl_ws1);
    gsl_integration_workspace_free(gsl_ws2);
    gsl_integration_workspace_free(gsl_ws3);

    return success;
}


real_t AvalancheSourceCH::distributionFunction(real_t xi, void *par){
    struct intFdistParams *params = (struct intFdistParams *) par;
    real_t p = params->pin;
    real_t mag = params->mag;
    real_t T = 50e4;
    real_t fac = DREAM::Constants::mc2inEV / (2 * M_PI * T);
    return mag * sqrt(fac*fac*fac) * exp(-p*p * DREAM::Constants::mc2inEV / (2 * T) - 5 * (1 - xi)); 
}

/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *AvalancheSourceCH::GetUnknownHandler(
    DREAM::FVM::Grid *g_fl, DREAM::FVM::Grid *g_hot, DREAM::FVM::Grid *g_re, 
    real_t n_re, const real_t n_tot, const real_t E_field
) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_RE, "0", g_re);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_F_HOT, "0", g_hot);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", g_fl);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_E_FIELD, "0", g_fl);

    real_t *temp_hot = new real_t[g_hot->GetNCells()];
    for (len_t ir = 0; ir < g_hot->GetNr(); ir++) {
        for (len_t j = 0; j < g_hot->GetNp2(ir); j++) {
            for (len_t i = 0; i < g_hot->GetNp1(ir); i++){
                real_t p = g_hot->GetMomentumGrid(ir)->GetP1(i);
                real_t xi = g_hot->GetMomentumGrid(ir)->GetP2(j);
                intFdistParams par_Fdist = {p, n_re};
                temp_hot[ir*g_hot->GetNp2(ir)*g_hot->GetNp1(ir) + j*g_hot->GetNp1(ir) + i] = distributionFunction(xi, &par_Fdist);
            }
        }
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_F_HOT, temp_hot); 

    real_t *temp_re  = new real_t[g_re->GetNCells()];
    for (len_t ir = 0; ir < g_re->GetNr(); ir++) {
        for (len_t j = 0; j < g_re->GetNp2(ir); j++) {
            for (len_t i = 0; i < g_re->GetNp1(ir); i++){
                real_t p = g_re->GetMomentumGrid(ir)->GetP1(i);
                real_t xi = g_re->GetMomentumGrid(ir)->GetP2(j);
                intFdistParams par_Fdist = {p, n_re};
                temp_re[ir*g_re->GetNp2(ir)*g_re->GetNp1(ir) + j*g_re->GetNp1(ir) + i] = distributionFunction(xi, &par_Fdist); 
            }
        }
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_F_RE, temp_re); 
    
    #define SETVAL_RE(NAME, v) do { \
            for (len_t i = 0; i < g_re->GetNCells(); i++) \
                temp_re[i] = (v); \
            uqh->SetInitialValue((NAME), temp_re); \
        } while (false)
    
    #define SETVAL_FL(NAME, v) do { \
            for (len_t i = 0; i < g_fl->GetNr(); i++) \
                temp_fl[i] = (v); \
            uqh->SetInitialValue((NAME), temp_fl); \
        } while (false)
    

    // Set electron quantities
    
    
    real_t *temp_fl = new real_t[g_fl->GetNr()];
    SETVAL_FL(DREAM::OptionConstants::UQTY_N_TOT, n_tot);
    SETVAL_FL(DREAM::OptionConstants::UQTY_E_FIELD, E_field);
    return uqh;
}


real_t AvalancheSourceCH::analyticIntegrandXi(real_t xi, void *par){
    struct intXiParams *params = (struct intXiParams *) par;
    real_t p = params->pin;
    real_t mag = params->mag;
    int QAG_KEY = params->qag_key;
    gsl_integration_workspace *gsl_ws3 = params->gslws3;
    real_t gamma = sqrt(p*p+1);
    real_t p_in = 2. * xi * p / (xi*xi * (gamma + 1.) - gamma + 1.);
    real_t gamma_in = (xi*xi * (gamma + 1.) + gamma - 1.) / (xi*xi * (gamma + 1.) - gamma + 1.);
    real_t gamma_in2 = gamma_in*gamma_in;
    real_t Sigma =  gamma_in2 / ((gamma_in2 - 1.) * (gamma - 1.)*(gamma - 1.) * (gamma_in - gamma)*(gamma_in - gamma)) 
                    * ((gamma_in - 1.)*(gamma_in - 1.) - (gamma - 1.) * (gamma_in - gamma) / (gamma_in*gamma_in) 
                        * (2. * gamma_in*gamma_in + 2. * gamma_in - 1. - (gamma - 1.) * (gamma_in - gamma)));

    
    real_t Fdist_neg, Fdist_pos, epsabs = 0, epsrel = 1e-12, lim = gsl_ws3->limit, error;
    gsl_function int_gsl_func;
    int_gsl_func.function = &(distributionFunction);
    intFdistParams par_Fdist = {p_in, mag};
    int_gsl_func.params = &par_Fdist;
    gsl_integration_qag(&int_gsl_func,-1,0,epsabs,epsrel,lim,QAG_KEY,gsl_ws3,&Fdist_neg, &error);
    gsl_integration_qag(&int_gsl_func, 0,1,epsabs,epsrel,lim,QAG_KEY,gsl_ws3,&Fdist_pos, &error);

    return p_in*p_in*p_in*p_in / xi * Sigma * (Fdist_neg + Fdist_pos) / 2;
}

real_t AvalancheSourceCH::analyticIntegrandP(real_t p, void *par){
    struct intPParams *params = (struct intPParams *) par;
    real_t p_max = params->pmax;
    real_t n_tot = params->ntot;
    real_t mag = params->mag;
    int QAG_KEY = params->qag_key;
    gsl_integration_workspace *gsl_ws2 = params->gslws2;
    gsl_integration_workspace *gsl_ws3 = params->gslws3;
    real_t gamma = sqrt(p*p + 1);
    real_t gamma_max = sqrt(p_max*p_max + 1);

    real_t ximin = sqrt((gamma - 1) / (gamma + 1) * (gamma_max + 1) / (gamma_max - 1));
    real_t ximax = sqrt(gamma / (gamma + 1));
    
    if (ximax <= ximin)
        return 0.;

    real_t intXi, epsabs = 0, epsrel = 1e-12, lim = gsl_ws2->limit, error;
    gsl_function int_gsl_func;
    int_gsl_func.function = &(analyticIntegrandXi);
    intXiParams par_Xi = {p, mag, QAG_KEY, gsl_ws3};
    int_gsl_func.params = &par_Xi;
    gsl_integration_qag(&int_gsl_func,ximin,ximax,epsabs,epsrel,lim, QAG_KEY,gsl_ws2,&intXi, &error);
    
    return 2 * M_PI * p*p * n_tot * (2 * M_PI * DREAM::Constants::r0 * DREAM::Constants::r0 * DREAM::Constants::c) * 1 / (p * gamma) * intXi;
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
    real_t successRelErrorThreshold = 0.01;

    len_t nr = 1;
    len_t np = 300;
    len_t nxi = 600;

    real_t pMin_hot = 0;
    real_t pMax_hot = 2;

    real_t pMin_re = 2;
    real_t pMax_re = 12;
    
    real_t B0 = 1.0;

    auto *grid_fl = InitializeFluidGrid(nr);
    auto *grid_hot = InitializeGridRCylPXi(nr,np,nxi,B0,pMin_hot,pMax_hot,true,pMax_re);
    auto *grid_re  = InitializeGridRCylPXi(nr,np,nxi,B0,pMin_re,pMax_re,true,pMax_re);

    real_t n_re = 1e15;
    real_t n_tot = 1e20;
    real_t E_field = 1;

    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(grid_fl, grid_hot, grid_re, n_re, n_tot, E_field);

    real_t pCutoff = 0.125;
    DREAM::AvalancheSourceCH *avaSourceTerm_hot = new DREAM::AvalancheSourceCH(
            grid_hot, uqh, pCutoff, 1.0, DREAM::AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, false, grid_re
        );
    
    DREAM::AvalancheSourceCH *avaSourceTerm_re = new DREAM::AvalancheSourceCH(
            grid_re, uqh, pCutoff, 1.0, DREAM::AvalancheSourceCH::CH_SOURCE_MODE_KINETIC, true
        );
    
    real_t *sourceVec_hot = new real_t[nr*np*nxi];
    real_t *sourceVec_re = new real_t[nr*np*nxi];
    
    real_t *deltas = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<np; i++){
            for(len_t j=0; j<nxi; j++){
                sourceVec_hot[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_hot->GetSourceFunction(ir,i,j);
                sourceVec_re[ir*np*nxi + j*np + i] = n_tot * avaSourceTerm_re->GetSourceFunction(ir,i,j);
            }
        }
        real_t sourceIntegralNumerical_hot = grid_hot->IntegralMomentumAtRadius(ir, sourceVec_hot);
        real_t sourceIntegralNumerical_re  =  grid_re->IntegralMomentumAtRadius(ir, sourceVec_re);

        real_t sourceIntegralAnalytic_hot, sourceIntegralAnalytic_re, epsabs = 0, epsrel = 1e-12, lim = gsl_ws1->limit, error;
        gsl_function int_gsl_func;
        int_gsl_func.function = &(analyticIntegrandP);
        intPParams par_P = {pMax_re, n_tot, n_re, QAG_KEY, gsl_ws2, gsl_ws3};
        int_gsl_func.params = &par_P;
        gsl_integration_qag(&int_gsl_func,pCutoff,pMax_hot,epsabs,epsrel,lim,QAG_KEY,gsl_ws1,&sourceIntegralAnalytic_hot, &error);
        gsl_integration_qag(&int_gsl_func,pMin_re,pMax_re,epsabs,epsrel,lim,QAG_KEY,gsl_ws1,&sourceIntegralAnalytic_re, &error);
        deltas[ir] = abs(sourceIntegralAnalytic_hot - sourceIntegralNumerical_hot) / sourceIntegralAnalytic_hot
                        + abs(sourceIntegralAnalytic_re - sourceIntegralNumerical_re) / sourceIntegralAnalytic_re;
        
    }


    bool success = true;
    for(len_t ir=0; ir<nr; ir++)
        if(deltas[ir]>successRelErrorThreshold){
            success = false;
            cout << "delta [rel error]: " << deltas[ir] << endl;
        }

    delete [] sourceVec_hot;
    delete [] sourceVec_re;
    return success;
}


