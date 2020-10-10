/**
 * Implementation of a radial grid in analytic toroidal magnetic geometry 
 * (with given magnetic-axis major radius, plasma minor radius and profiles of elongation, 
 *  triangularity and Shafranov shift as well as reference poloidal flux profile)
 */

#include <cmath>
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>

using namespace DREAM::FVM;

/**
 * Constructor.
 * nx: Number of radial grid points.
 * G: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 * x0: Value of inner radial flux grid point.
 * xa: Value of outer radial flux grid point.
 */

AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const len_t nr,  real_t r0,  real_t ra, real_t R0, len_t ntheta_interp,
    real_t *rProfiles, len_t nrProfiles, real_t *Gs, real_t *psi_p0s,
             real_t *kappas, real_t *deltas, real_t *Deltas
) : RadialGridGenerator(nr), rMin(r0), rMax(ra) {
    this->R0             = R0;
    this->ntheta_interp  = ntheta_interp;
    this->GsProvided     = Gs;
    this->psisProvided   = psi_p0s;
    this->kappasProvided = kappas;
    // note: deltas and Deltas should vanish at rProfiles = 0 (physical constraint, fine numerically regardless)
    this->deltasProvided = deltas;
    this->DeltasProvided = Deltas;
    this->nrProfiles     = nrProfiles;
    this->rProfilesProvided = rProfiles;

    isUpDownSymmetric = true;
    spline_x = gsl_spline_alloc(gsl_interp_steffen, nrProfiles);
    gsl_acc  = gsl_interp_accel_alloc();
}

AnalyticBRadialGridGenerator::~AnalyticBRadialGridGenerator(){
    gsl_spline_free (spline_x);
    gsl_interp_accel_free (gsl_acc);
    DeallocateShapeProfiles();
}

/**
 * (Re-)builds the given radial grid.
 *
 * rGrid: Radial grid to re-build.
 */
bool AnalyticBRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    r    = new real_t[GetNr()];
    r_f  = new real_t[GetNr()+1];
    
    real_t
        *dr   = new real_t[GetNr()],
        *dr_f = new real_t[GetNr()-1];

    // Construct flux grid
    for (len_t i = 0; i < GetNr(); i++)
        dr[i] = (rMax - rMin) / GetNr();

    for (len_t i = 0; i < GetNr()+1; i++)
        r_f[i] = rMin + i*dr[0];

    // Construct cell grid
    for (len_t i = 0; i < GetNr(); i++)
        r[i] = 0.5 * (r_f[i+1] + r_f[i]);

    for (len_t i = 0; i < GetNr()-1; i++)
        dr_f[i] = r[i+1] - r[i];

    rGrid->Initialize(r, r_f, dr, dr_f);

    DeallocateShapeProfiles();
    InterpolateInputProfileToGrid(r,r_f,BtorGOverR0,GPrime,BtorGOverR0_f,GPrime_f,GsProvided);
    InterpolateInputProfileToGrid(r,r_f,psi,psiPrimeRef,psi_f,psiPrimeRef_f,psisProvided);
    InterpolateInputProfileToGrid(r,r_f,kappa,kappaPrime,kappa_f,kappaPrime_f,kappasProvided);
    InterpolateInputProfileToGrid(r,r_f,delta,deltaPrime,delta_f,deltaPrime_f,deltasProvided);
    InterpolateInputProfileToGrid(r,r_f,Delta,DeltaPrime,Delta_f,DeltaPrime_f,DeltasProvided);
    rGrid->SetReferenceMagneticFieldData(
        BtorGOverR0, BtorGOverR0_f, psiPrimeRef, psiPrimeRef_f, R0
    );

    this->isBuilt = true;
    return true;
}

/**
 *  Numerically differentiates the input (lambda) function F(r) with respect to its argument
 */
real_t AnalyticBRadialGridGenerator::diffFunc(real_t r, std::function<real_t(real_t)> F){
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h = sqrteps * ( 1 + abs(r) ); 
    return (F(r+h/2)-F(r-h/2))/h;
}

/**
 * Evaluates the local major radius at radial grid point ir and poloidal angle theta 
 */
real_t AnalyticBRadialGridGenerator::ROverR0AtTheta(const len_t ir, const real_t theta) {
    if(isinf(R0))
        return 1;
    else
        return 1 + (Delta[ir] + r[ir]*cos(theta + delta[ir]*sin(theta)))/R0;
}
// Same as ROverR0AtTheta but evaluated on the radial flux grid
real_t AnalyticBRadialGridGenerator::ROverR0AtTheta_f(const len_t ir, const real_t theta) {
    if(isinf(R0))
        return 1;
    else
        return 1 + (Delta_f[ir] + r_f[ir]*cos(theta + delta_f[ir]*sin(theta)))/R0;
}


// Evaluates the spatial Jacobian normalized to r*R
real_t AnalyticBRadialGridGenerator::normalizedJacobian(const len_t ir, const real_t theta){
    real_t ct = cos(theta);
    real_t st = sin(theta);
    return ( kappa[ir]*cos(delta[ir]*st) + kappa[ir]*DeltaPrime[ir]*ct
        + st*sin(theta+delta[ir]*st) * ( r[ir]*kappaPrime[ir] +
        ct * (  delta[ir]*kappa[ir] + r[ir]* delta[ir]*kappaPrime[ir] - r[ir]*kappa[ir]*deltaPrime[ir] ) ) );
}
// Evaluates the spatial Jacobian normalized to r*R on the radial flux grid
real_t AnalyticBRadialGridGenerator::normalizedJacobian_f(const len_t ir, const real_t theta){
    real_t ct = cos(theta);
    real_t st = sin(theta);
    return ( kappa_f[ir]*cos(delta_f[ir]*st) + kappa_f[ir]*DeltaPrime_f[ir]*ct
        + st*sin(theta+delta_f[ir]*st) * ( r_f[ir]*kappaPrime_f[ir] +
        ct * (  delta_f[ir]*kappa_f[ir] + r_f[ir]* delta_f[ir]*kappaPrime_f[ir] - r_f[ir]*kappa_f[ir]*deltaPrime_f[ir] ) ) );
}

/**
 * Evaluates the spatial Jacobian normalized to R0 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::JacobianAtTheta(const len_t ir, const real_t theta){
    return r[ir]*ROverR0AtTheta(ir,theta) * normalizedJacobian(ir,theta);
}
// Same as JacobianAtTheta but evaluated on the radial flux grid
real_t AnalyticBRadialGridGenerator::JacobianAtTheta_f(const len_t ir, const real_t theta){
    return r_f[ir]*ROverR0AtTheta_f(ir,theta) * normalizedJacobian_f(ir,theta);
}

/**
 * Evaluates |nabla r|^2 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::NablaR2AtTheta(const len_t ir, const real_t theta){
    real_t st = sin(theta);
    real_t ct = cos(theta);
    real_t JOverRr = normalizedJacobian(ir,theta);
    return (kappa[ir]*kappa[ir] * ct*ct + (1+delta[ir]*ct)*(1+delta[ir]*ct) 
                * sin(theta+delta[ir]*st)*sin(theta+delta[ir]*st) )  / (JOverRr*JOverRr); 
}
/**
 * Evaluates |nabla r|^2 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::NablaR2AtTheta_f(const len_t ir, const real_t theta){
    real_t st = sin(theta);
    real_t ct = cos(theta);
    real_t JOverRr = normalizedJacobian_f(ir,theta);
    return (kappa_f[ir]*kappa_f[ir] * ct*ct + (1+delta_f[ir]*ct)*(1+delta_f[ir]*ct) 
                * sin(theta+delta_f[ir]*st)*sin(theta+delta_f[ir]*st) )  / (JOverRr*JOverRr); 
}

/**
 * Interpolates input shape-parameter profiles (kappa, delta, ...) which are defined on 
 * input rProfilesProvided array to the r and r_f grids
 */
void AnalyticBRadialGridGenerator::InterpolateInputProfileToGrid(real_t *r, real_t *r_f, real_t *&x,real_t *&xPrime, real_t *&x_f, real_t *&xPrime_f,real_t *xProvided){
    x        = new real_t[GetNr()];
    xPrime   = new real_t[GetNr()];
    x_f      = new real_t[GetNr()+1];
    xPrime_f = new real_t[GetNr()+1];

    gsl_spline_init(spline_x, rProfilesProvided, xProvided, nrProfiles);
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h;
    for (len_t ir=0; ir<GetNr(); ir++){
        x[ir]      = gsl_spline_eval(spline_x,r[ir],gsl_acc);
        
        h = sqrteps * ( 1 + abs(r[ir]) );
        if (ir==GetNr()-1)
            h = -h;
        xPrime[ir] = (gsl_spline_eval(spline_x,r[ir]+h,gsl_acc) - x[ir])/h;
    }
    for (len_t ir=0; ir<GetNr()+1; ir++){
        x_f[ir]      = gsl_spline_eval(spline_x,r_f[ir],gsl_acc);
        
        h = sqrteps * ( 1 + abs(r_f[ir]) );
        if (ir==GetNr())
            h = -h;
        xPrime_f[ir] = (gsl_spline_eval(spline_x,r_f[ir]+h,gsl_acc) - x_f[ir])/h;
    }
}

/**
 * Deallocates shape profiles that have been interpolated to the grid.
 */
void AnalyticBRadialGridGenerator::DeallocateShapeProfiles(){
    if (psi==nullptr)
        return;

    delete [] psi;
    delete [] kappa;
    delete [] delta;
    delete [] Delta;
    delete [] psi_f;
    delete [] kappa_f;
    delete [] delta_f;
    delete [] Delta_f;
    delete [] GPrime;
    delete [] kappaPrime;
    delete [] deltaPrime;
    delete [] DeltaPrime;
    delete [] GPrime_f;
    delete [] kappaPrime_f;
    delete [] deltaPrime_f;
    delete [] DeltaPrime_f;
}

