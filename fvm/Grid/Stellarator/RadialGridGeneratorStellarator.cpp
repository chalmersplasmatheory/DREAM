
#include <algorithm>
#include <vector>
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/Stellarator/RadialGridGeneratorStellarator.hpp"
#include <gsl/gsl_errno.h>

using namespace std;
using namespace DREAM::FVM;

RadialGridGeneratorStellarator::RadialGridGeneratorStellarator(const len_t nr) 
        : RadialGridGenerator(nr) {
    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2; // TODO: Ok?
    gsl_multi_fmin = gsl_multimin_fminimizer_alloc(T, 2);
}


RadialGridGeneratorStellarator::~RadialGridGeneratorStellarator(){
    gsl_multimin_fminimizer_free(gsl_multi_fmin);
}


/** TODO: Is this needed?
 * Rebuilds magnetic field data and stores all quantities needed for flux surface and bounce averages.
 */
void RadialGridGeneratorStellarator::RebuildJacobians(RadialGridStellarator *rGrid) {
    nr = rGrid->GetNr();

    Bmin   = new real_t[GetNr()];
    Bmax   = new real_t[GetNr()];
    Bmin_f = new real_t[GetNr()+1];
    Bmax_f = new real_t[GetNr()+1];
    theta_Bmin     = new real_t[GetNr()];
    theta_Bmax     = new real_t[GetNr()];
    theta_Bmin_f   = new real_t[GetNr()+1];
    theta_Bmax_f   = new real_t[GetNr()+1];
    phi_Bmin     = new real_t[GetNr()];
    phi_Bmax     = new real_t[GetNr()];
    phi_Bmin_f   = new real_t[GetNr()+1];
    phi_Bmax_f   = new real_t[GetNr()+1];
    xi0TrappedBoundary   = new real_t[GetNr()];
    xi0TrappedBoundary_f = new real_t[GetNr()+1];
    std::array<real_t,2> extremum;
    for (len_t ir = 0; ir < GetNr(); ir++){
        extremum = getThetaPhi_Bmin(ir);
        theta_Bmin[ir] = extremum[0];
        phi_Bmin[ir] = extremum[1];

        extremum = getThetaPhi_Bmax(ir);
        theta_Bmax[ir] = extremum[0];
        phi_Bmax[ir] = extremum[1];
        
        Bmin[ir] = BAtThetaPhi(ir,theta_Bmin[ir], phi_Bmin[ir]);
        Bmax[ir] = BAtThetaPhi(ir,theta_Bmax[ir], phi_Bmax[ir]);
        xi0TrappedBoundary[ir] = 0; // Used for BA in E_ceff, simplified => no trapped particles
        
    }
    for (len_t ir = 0; ir < GetNr()+1; ir++){
        extremum = getThetaPhi_Bmin_f(ir);
        theta_Bmin_f[ir] = extremum[0];
        phi_Bmin_f[ir] = extremum[1];

        extremum = getThetaPhi_Bmax_f(ir);
        theta_Bmax_f[ir] = extremum[0];
        phi_Bmax_f[ir] = extremum[1];

        Bmin_f[ir] = BAtThetaPhi_f(ir,theta_Bmin_f[ir], phi_Bmin_f[ir]);
        Bmax_f[ir] = BAtThetaPhi_f(ir,theta_Bmax_f[ir], phi_Bmax_f[ir]);
        xi0TrappedBoundary_f[ir] = 0;  // Used for BA in E_ceff, simplified => no trapped particles
    }
    rGrid->SetMagneticExtremumData(
        Bmin, Bmin_f, Bmax, Bmax_f, 
        theta_Bmin, theta_Bmin_f, theta_Bmax, theta_Bmax_f,
        xi0TrappedBoundary, xi0TrappedBoundary_f
    );
    delete [] phi_Bmin;
    delete [] phi_Bmax;
}

// The remaining functions are related to determining theta_Bmin and theta_Bmax
// with a gsl fmin algorithm
struct EvalBParams {len_t ir; RadialGridGeneratorStellarator* rgg; int_t sgn;};
real_t gslEvalB(const gsl_vector *v, void *par){
    real_t theta = gsl_vector_get(v, 0);
    real_t phi = gsl_vector_get(v, 1);
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtThetaPhi(params->ir,theta,phi);
}
real_t gslEvalB_f(const gsl_vector *v, void *par){
    real_t theta = gsl_vector_get(v, 0);
    real_t phi = gsl_vector_get(v, 1);
    EvalBParams *params = (EvalBParams *) par;
    return params->sgn*params->rgg->BAtThetaPhi_f(params->ir,theta, phi);
}


// Tolerances and max number of iterations for gsl magnetic field minimizer
const len_t MAX_NUM_ITER = 30;
const real_t EPSABS = 1e-6;
const real_t EPSREL = 0;
const real_t STEP = 2*M_PI / 100; 
/**
 * Finds the extremum of the magnetic field on the interval [0,2*pi]. 
 * If sgn=1, returns the minimum.
 * If sgn=-1, returns the maximum.
 */
std::array<real_t,2> RadialGridGeneratorStellarator::FindMagneticFieldExtremumStellarator(
    len_t ir, int_t sgn, fluxGridType fluxGridType
) { 
    real_t theta_guess = 0, phi_guess = 0;
    real_t B_opt = sgn*std::numeric_limits<real_t>::infinity();
    real_t B;

    real_t phi_max = 2*M_PI;
    if (this->nfp > 0)
        phi_max = M_PI / this->nfp;
    
    for (real_t theta=0; theta<2*M_PI; theta+=STEP){
        for (real_t phi=0; phi<phi_max; phi+=STEP){
            if (fluxGridType == FLUXGRIDTYPE_DISTRIBUTION) {
                B = BAtThetaPhi(ir, theta, phi);
            } else {
                B = BAtThetaPhi_f(ir, theta, phi);
            }
            if (sgn*B < sgn*B_opt){
                B_opt = B;
                theta_guess = theta;
                phi_guess = phi;
            }
        }
    }

    gsl_vector *guess = gsl_vector_alloc(2);
    gsl_vector_set(guess, 0, theta_guess);
    gsl_vector_set(guess, 1, phi_guess);

    gsl_vector *step = gsl_vector_alloc(2);
    gsl_vector_set_all(step, STEP / 2); // TODO: OK Step size?
	
    EvalBParams params = {ir, this, sgn};
    gsl_multimin_function gsl_func;
    gsl_func.n = 2;
    if(fluxGridType == FLUXGRIDTYPE_DISTRIBUTION)
        gsl_func.f = &(gslEvalB); 
    else
        gsl_func.f = &(gslEvalB_f); 
    gsl_func.params = &(params); // TODO: Should this be without &?

    
    // otherwise, find extremum with fmin algorithm
    gsl_multimin_fminimizer_set( 
        gsl_multi_fmin, &gsl_func, guess, step 
    );

    int status;
    real_t size;
    for(len_t iter=0; iter<MAX_NUM_ITER; iter++){
        gsl_multimin_fminimizer_iterate(gsl_multi_fmin);

        size = gsl_multimin_fminimizer_size(gsl_multi_fmin);
        status = gsl_multimin_test_size(size, EPSABS);
        if(status == GSL_SUCCESS)
            break;
    }

    gsl_vector_free(guess);
    gsl_vector_free(step);

    real_t theta = gsl_vector_get(gsl_multi_fmin->x, 0);
    if(theta < 2*EPSABS || theta > (2*M_PI - 2*EPSABS))
        theta = 0;
    else if (fabs(M_PI-theta) < 2*EPSABS)
        theta = M_PI;
    real_t phi   = gsl_vector_get(gsl_multi_fmin->x, 1);
    
    std::array<real_t,2> extremum;
    extremum[0] = theta;
    extremum[1] = phi;
    return extremum;
}
