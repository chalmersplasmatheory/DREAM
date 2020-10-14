
#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <gsl/gsl_errno.h>

using namespace std;
using namespace DREAM::FVM;

RadialGridGenerator::RadialGridGenerator(const len_t nr) : nr(nr) {
    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_fmin = gsl_min_fminimizer_alloc(T); 

}


RadialGridGenerator::~RadialGridGenerator(){
    gsl_min_fminimizer_free(gsl_fmin);
}


/**
 * Rebuilds magnetic field data and stores all quantities needed for flux surface and bounce averages.
 */
void RadialGridGenerator::RebuildJacobians(RadialGrid *rGrid) {
    nr = rGrid->GetNr();

    Bmin   = new real_t[GetNr()];
    Bmax   = new real_t[GetNr()];
    Bmin_f = new real_t[GetNr()+1];
    Bmax_f = new real_t[GetNr()+1];
    theta_Bmin     = new real_t[GetNr()];
    theta_Bmax     = new real_t[GetNr()];
    theta_Bmin_f   = new real_t[GetNr()+1];
    theta_Bmax_f   = new real_t[GetNr()+1];
    for (len_t ir = 0; ir < GetNr(); ir++){
        theta_Bmin[ir] = getTheta_Bmin(ir);
        theta_Bmax[ir] = getTheta_Bmax(ir);
        Bmin[ir] = BAtTheta(ir,theta_Bmin[ir]);
        Bmax[ir] = BAtTheta(ir,theta_Bmax[ir]);        
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        theta_Bmin_f[ir] = getTheta_Bmin_f(ir);
        theta_Bmax_f[ir] = getTheta_Bmax_f(ir);
        Bmin_f[ir] = BAtTheta_f(ir,theta_Bmin_f[ir]);
        Bmax_f[ir] = BAtTheta_f(ir,theta_Bmax_f[ir]);
    }
    rGrid->SetMagneticExtremumData(
        Bmin, Bmin_f, Bmax, Bmax_f, 
        theta_Bmin, theta_Bmin_f, theta_Bmax, theta_Bmax_f
    );
}


// Evaluates the magnetic field strength at radial index ir 
// on the distribution grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta(const len_t ir, const real_t theta) {
    real_t ct = cos(theta);
    real_t st = sin(theta);
    return BAtTheta(ir,theta,ct,st);
}
// Evaluates the magnetic field strength at radial index ir 
// on the distribution grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta(const len_t ir, const real_t theta, const real_t ct, const real_t st) {
    real_t ROverR0 = ROverR0AtTheta(ir,theta,ct,st);
    real_t Btor = BtorGOverR0[ir]/ROverR0;
    real_t Bpol = 0;
    if(psiPrimeRef[ir])
        Bpol = sqrt(NablaR2AtTheta(ir,theta,ct,st))*psiPrimeRef[ir]/ROverR0;  
    return sqrt(Btor*Btor+Bpol*Bpol);
}
// Evaluates the magnetic field strength at radial index ir 
// on the radial flux grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta_f(const len_t ir, const real_t theta) {
    real_t ct = cos(theta);
    real_t st = sin(theta);
    return BAtTheta_f(ir,theta,ct,st);
}
// Evaluates the magnetic field strength at radial index ir 
// on the radial flux grid and poloidal angle theta
real_t RadialGridGenerator::BAtTheta_f(const len_t ir, const real_t theta, const real_t ct, const real_t st) {
    real_t ROverR0 = ROverR0AtTheta_f(ir,theta,ct,st);
    real_t Btor = BtorGOverR0_f[ir]/ROverR0;
    real_t Bpol = 0;
    if(psiPrimeRef_f[ir])
        Bpol = sqrt(NablaR2AtTheta_f(ir,theta,ct,st))*psiPrimeRef_f[ir]/ROverR0;  
    return sqrt(Btor*Btor+Bpol*Bpol);
}





// The remaining functions are related to determining theta_Bmin and theta_Bmax
// with a gsl fmin algorithm

struct EvalBParams {len_t ir; RadialGridGenerator* rgg; int_t sgn;};
real_t gslEvalB(real_t theta, void *par){
    EvalBParams *params = (EvalBParams *) par;
    len_t ir = params->ir;
    int_t sgn = params->sgn;
    RadialGridGenerator *rgg = params->rgg;
    return sgn*rgg->BAtTheta(ir,theta);
}
real_t gslEvalB_f(real_t theta, void *par){
    EvalBParams *params = (EvalBParams *) par;
    len_t ir = params->ir;
    int_t sgn = params->sgn;
    RadialGridGenerator *rgg = params->rgg;
    return sgn*rgg->BAtTheta_f(ir,theta);
}


// Tolerances and max number of iterations for gsl magnetic field minimizer
const len_t MAX_NUM_ITER = 30;
const real_t EPSABS = 1e-6;
const real_t EPSREL = 0;
/**
 * Finds the extremum of the magnetic field on the interval [0,pi]. 
 * If sgn=1, returns the minimum.
 * If sgn=-1, returns the maximum.
 */
real_t RadialGridGenerator::FindMagneticFieldExtremum(len_t ir, int_t sgn, fluxGridType fluxGridType){
    if(!isUpDownSymmetric)
        throw FVMException(
            "RadialGridGenerator: for magnetic geometry which is not up-down symmetric,"
            " must override getTheta functions."
        );
    EvalBParams params = {ir, this, sgn};
    gsl_function gsl_func;
    if(fluxGridType == FLUXGRIDTYPE_RADIAL)
        gsl_func.function = &(gslEvalB_f);
    else
        gsl_func.function = &(gslEvalB);
    gsl_func.params = &(params);

    real_t theta_lim_lower = 0;
    real_t theta_lim_upper = M_PI;
    real_t theta_guess;
    if(sgn==1){
        // if B has a local minimum in theta=0, return 0
        theta_guess = 10*EPSABS;
        if(gsl_func.function(theta_guess,gsl_func.params) >= gsl_func.function(theta_lim_lower,gsl_func.params))
            return 0;
    } else { 
        // if B has a local maximum in theta=pi, return pi
        theta_guess=M_PI-10*EPSABS;
        if(gsl_func.function(theta_guess,gsl_func.params) >= gsl_func.function(theta_lim_upper,gsl_func.params))
            return M_PI;
    }
    // otherwise, find extremum with fmin algorithm
    gsl_min_fminimizer_set(
        gsl_fmin, &gsl_func, theta_guess, theta_lim_lower, theta_lim_upper
    );

    int status;
    for(len_t iter=0; iter<MAX_NUM_ITER; iter++){
        gsl_min_fminimizer_iterate (gsl_fmin);
        real_t x_lo = gsl_min_fminimizer_x_lower (gsl_fmin);
        real_t x_up = gsl_min_fminimizer_x_upper (gsl_fmin);
        status = gsl_min_test_interval(x_lo, x_up, EPSABS, EPSREL);
        if(status == GSL_SUCCESS)
            break;
    }
    real_t extremum = gsl_min_fminimizer_x_minimum(gsl_fmin); 
    if(extremum < 2*EPSABS)
        return 0;
    else if (abs(M_PI-extremum) < 2*EPSABS)
        return M_PI;
    else
        return extremum; 
}
