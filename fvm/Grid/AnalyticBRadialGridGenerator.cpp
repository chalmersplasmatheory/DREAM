/**
 * Implementation of a radial grid in analytic toroidal magnetic geometry 
 * (with given magnetic-axis major radius, plasma minor radius and profiles of elongation, 
 *  triangularity and Shafranov shift as well as reference poloidal flux profile)
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"

#include <math.h>

using namespace DREAM::FVM;
using namespace std;

/**
 * Constructor.
 * nr: Number of radial grid points.
 * G: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 * r0: Value of inner radial flux grid point.
 * ra: Value of outer radial flux grid point.
 */
AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const len_t nr,  real_t r0,  real_t ra, real_t R0, len_t ntheta_interp,
     struct shape_profiles *profiles
) : RadialGridGenerator(nr), rMin(r0), rMax(ra), providedProfiles(profiles) {
    this->R0             = R0;
    this->ntheta_interp  = ntheta_interp;

    R0IsInf = isinf(R0);
    constructSplines(profiles);
    isUpDownSymmetric = true;
}

/**
 * Constructor.
 * 
 * r_f_input: Grid points on the radial flux grid (e.g. the cell edges)
 * nr: Number of radial grid points.
 * G: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 */
AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const real_t *r_f_input, const len_t nr, real_t R0, len_t ntheta_interp,
     struct shape_profiles *profiles
) : RadialGridGenerator(nr), rMin(r_f_input[0]), rMax(r_f_input[nr]), providedProfiles(profiles) {
    this->R0             = R0;
    this->ntheta_interp  = ntheta_interp;

    this->rf_provided = new real_t[nr+1];
    for(len_t i=0; i<nr+1; i++)
        this->rf_provided[i] = r_f_input[i];

    delete [] r_f_input;

    constructSplines(profiles);
    isUpDownSymmetric = true;
}



/**
 * Destructor
 */
AnalyticBRadialGridGenerator::~AnalyticBRadialGridGenerator(){
    if (this->spline_G != nullptr) {
        gsl_spline_free(spline_G);
        gsl_interp_accel_free(gsl_acc_G);
    }
    if (this->spline_psi != nullptr) {
        gsl_spline_free(spline_psi);
        gsl_interp_accel_free(gsl_acc_psi);
    }
    if (this->spline_kappa != nullptr) {
        gsl_spline_free(spline_kappa);
        gsl_interp_accel_free(gsl_acc_kappa);
    }
    if (this->spline_delta != nullptr) {
        gsl_spline_free(spline_delta);
        gsl_interp_accel_free(gsl_acc_delta);
    }
    if (this->spline_Delta != nullptr) {
        gsl_spline_free(spline_Delta);
        gsl_interp_accel_free(gsl_acc_Delta);
    }

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

    // if r_f has been provided to constructor, set specified grid, otherwise uniform
    if(rf_provided==nullptr){
        for (len_t i = 0; i < GetNr(); i++)
            dr[i] = (rMax - rMin) / GetNr();
        for (len_t i = 0; i < GetNr()+1; i++)
            r_f[i] = rMin + i*dr[0];
    } else { 
        for(len_t i=0; i < GetNr()+1; i++)
            r_f[i] = rf_provided[i];
        for (len_t i = 0; i < GetNr(); i++)
            dr[i] = r_f[i+1] - r_f[i];
    }
    delete [] rf_provided;
    rf_provided = nullptr; // upon next rebuild, create uniform grid at the new resolution

    // Construct cell grid
    for (len_t i = 0; i < GetNr(); i++)
        r[i] = 0.5 * (r_f[i+1] + r_f[i]);

    for (len_t i = 0; i < GetNr()-1; i++)
        dr_f[i] = r[i+1] - r[i];

    rGrid->Initialize(r, r_f, dr, dr_f);

    struct shape_profiles *pp = this->providedProfiles;

    DeallocateShapeProfiles();
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nG,     pp->G,     spline_G,     gsl_acc_G,     &BtorGOverR0, &GPrime,      &BtorGOverR0_f, &GPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->npsi,   pp->psi,   spline_psi,   gsl_acc_psi,   &psi,         &psiPrimeRef, &psi_f,         &psiPrimeRef_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nkappa, pp->kappa, spline_kappa, gsl_acc_kappa, &kappa,       &kappaPrime,  &kappa_f,       &kappaPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->ndelta, pp->delta, spline_delta, gsl_acc_delta, &delta,       &deltaPrime,  &delta_f,       &deltaPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nDelta, pp->Delta, spline_Delta, gsl_acc_Delta, &Delta,       &DeltaPrime,  &Delta_f,       &DeltaPrime_f);
    if(r_f[0]==0) // standard situation
        psiPrimeRef_f[0] = 0; // no poloidal field at r=0 since no toroidal current is enclosed
    rGrid->SetReferenceMagneticFieldData(
        BtorGOverR0, BtorGOverR0_f, psiPrimeRef, psiPrimeRef_f, R0
    );

    this->isBuilt = true;
    return true;
}


/**
 * Evaluates the local major radius at radial grid point ir and poloidal angle theta 
 */
real_t AnalyticBRadialGridGenerator::ROverR0AtTheta(const len_t ir, const real_t theta) {
    if(R0IsInf)
        return 1;
    else
        return 1 + (Delta[ir] + r[ir]*cos(theta + delta[ir]*sin(theta)))/R0;
}
// Same as ROverR0AtTheta but evaluated on the radial flux grid
real_t AnalyticBRadialGridGenerator::ROverR0AtTheta_f(const len_t ir, const real_t theta) {
    if(R0IsInf || (r_f[ir]==0 && Delta_f[ir]==0) )
        return 1;
    else
        return 1 + (Delta_f[ir] + r_f[ir]*cos(theta + delta_f[ir]*sin(theta)))/R0;
}


// Evaluates the spatial Jacobian normalized to r*R
real_t AnalyticBRadialGridGenerator::normalizedJacobian(const len_t ir, const real_t theta, real_t ct, real_t st){
    return kappa[ir]*cos(delta[ir]*st) + kappa[ir]*DeltaPrime[ir]*ct
        + st*sin(theta+delta[ir]*st) * ( r[ir]*kappaPrime[ir] +
        ct * (  delta[ir]*kappa[ir] + r[ir]* delta[ir]*kappaPrime[ir]
               - r[ir]*kappa[ir]*deltaPrime[ir] ) ) ;
}
// Evaluates the spatial Jacobian normalized to r*R on the radial flux grid
real_t AnalyticBRadialGridGenerator::normalizedJacobian_f(const len_t ir, const real_t theta, real_t ct, real_t st){
    return kappa_f[ir]*cos(delta_f[ir]*st) + kappa_f[ir]*DeltaPrime_f[ir]*ct
        + st*sin(theta+delta_f[ir]*st) * ( r_f[ir]*kappaPrime_f[ir] +
        ct * (  delta_f[ir]*kappa_f[ir] + r_f[ir]* delta_f[ir]*kappaPrime_f[ir]
               - r_f[ir]*kappa_f[ir]*deltaPrime_f[ir] ) ) ;
}

/**
 * Evaluates the spatial Jacobian normalized to R0 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::JacobianAtTheta(const len_t ir, const real_t theta){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);

    return r[ir]*ROverR0AtTheta(ir,theta) * normalizedJacobian(ir,theta,ct,st);
}
// Same as JacobianAtTheta but evaluated on the radial flux grid
real_t AnalyticBRadialGridGenerator::JacobianAtTheta_f(const len_t ir, const real_t theta){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);

    return r_f[ir]*ROverR0AtTheta_f(ir,theta) * normalizedJacobian_f(ir,theta,ct,st);
}

/**
 * Evaluates |nabla r|^2 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::NablaR2AtTheta(const len_t ir, const real_t theta){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);

    real_t sdt = sin(theta+delta[ir]*st);
    real_t cdt = 1+delta[ir]*ct;
    real_t JOverRr = normalizedJacobian(ir,theta,ct,st);
    return (kappa[ir]*kappa[ir] * ct*ct + cdt * cdt
                * sdt * sdt ) / (JOverRr*JOverRr);
}
/**
 * Evaluates |nabla r|^2 at radial grid point ir and poloidal angle theta on the radial flux grid
 */
real_t AnalyticBRadialGridGenerator::NablaR2AtTheta_f(const len_t ir, const real_t theta){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);

    real_t sdt = sin(theta+delta_f[ir]*st);
    real_t cdt = 1+delta_f[ir]*ct;
    real_t JOverRr = normalizedJacobian_f(ir,theta,ct,st);
    return (kappa_f[ir]*kappa_f[ir] * ct * ct + cdt * cdt 
                * sdt * sdt)  / (JOverRr*JOverRr); 
}

/**
 * The method evaluates all geometric quantities at radial grid point ir and poloidal angle theta
 */
void AnalyticBRadialGridGenerator::EvaluateGeometricQuantities(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);
    real_t sdt = 0.0;
    real_t cdt = 1.0;
    if(delta[ir])
        sincos(delta[ir]*st,&sdt,&cdt);

    real_t stdt = st*cdt+sdt*ct; // = sin(theta + delta*sin(theta))
    real_t ctdt = ct*cdt-st*sdt; // = cos(theta + delta*sin(theta))

    real_t rk = r[ir]*kappaPrime[ir];
    real_t JOverRr = kappa[ir]*(cdt + DeltaPrime[ir]*ct)
        + st*stdt * ( rk + ct * (  delta[ir]*(kappa[ir] + rk)
               - r[ir]*kappa[ir]*deltaPrime[ir] ) ) ;
    
    ROverR0 = 1;
    if(!R0IsInf)
        ROverR0 += (Delta[ir] + r[ir]*ctdt)/R0;
    
    Jacobian = r[ir] * ROverR0 * JOverRr;

    real_t deltaTerm = (1.0+delta[ir]*ct)*stdt;
    real_t kappaTerm = kappa[ir]*ct;
    NablaR2 = (kappaTerm*kappaTerm + deltaTerm*deltaTerm) / (JOverRr*JOverRr);
    
    real_t Btor = BtorGOverR0[ir]/ROverR0;
    real_t BpolSq = 0;
    if(psiPrimeRef[ir]){
        real_t Bp0 = psiPrimeRef[ir]/(2*M_PI*ROverR0);
        BpolSq = NablaR2*Bp0*Bp0;  
    }
    B = sqrt(Btor*Btor+BpolSq);
}

/**
 * Same as EvaluateGeometricQuantities, but on the radial flux grid
 */
void AnalyticBRadialGridGenerator::EvaluateGeometricQuantities_fr(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2){
    real_t ct = 0; // cos(theta)
    real_t st = 0; // sin(theta)
    sincos(theta,&st,&ct);
    real_t sdt = 0.0;
    real_t cdt = 1.0;
    if(delta_f[ir])
        sincos(delta_f[ir]*st,&sdt,&cdt);

    real_t stdt = st*cdt+sdt*ct;  // = sin(theta + delta*sin(theta))
    real_t ctdt = ct*cdt-st*sdt;  // = cos(theta + delta*sin(theta))

    real_t rk = r_f[ir]*kappaPrime_f[ir];
    real_t JOverRr = kappa_f[ir]*(cdt + DeltaPrime_f[ir]*ct)
        + st*stdt * ( rk + ct * (  delta_f[ir]*(kappa_f[ir] + rk)
               - r_f[ir]*kappa_f[ir]*deltaPrime_f[ir] ) ) ;

    ROverR0 = 1;
    if(!R0IsInf)
        ROverR0 += (Delta_f[ir] + r_f[ir]*ctdt)/R0;

    Jacobian = r_f[ir] * ROverR0 * JOverRr;

    real_t deltaTerm = (1+delta_f[ir]*ct)*stdt;
    real_t kappaTerm = kappa_f[ir] * ct;
    NablaR2 = (kappaTerm*kappaTerm +  deltaTerm*deltaTerm)  / (JOverRr*JOverRr);
    
    real_t Btor = BtorGOverR0_f[ir]/ROverR0;
    real_t BpolSq = 0;
    if(psiPrimeRef_f[ir]){
        real_t Bp0 = psiPrimeRef_f[ir]/(2*M_PI*ROverR0);
        BpolSq = NablaR2*Bp0*Bp0;  
    }
    B = sqrt(Btor*Btor+BpolSq);
}


/**
 * Interpolates input shape-parameter profiles (kappa, delta, ...) which are defined on 
 * input rProfilesProvided array to the r and r_f grids
 *
 * nr:        Number of grid points on output grid.
 * r:         Output radial grid.
 * r_f:       Output radial flux grid.
 * nProvided: Number of grid points in which the shape parameter is specified.
 * rProvided: Radial grid on which the shape parameter is provided.
 * xProvided: Provided shape parameter.
 *
 * OUTPUT
 * x:         Shape parameter evaluated on 'r'.
 * xPrime:
 * x_f:       Shape parameter evaluated on 'r_f'.
 * xPrime_f:
 */
void AnalyticBRadialGridGenerator::InterpolateInputProfileToGrid(
    const len_t nr, const real_t *r, const real_t *r_f,
    const len_t nProvided, const real_t *xProvided,
    gsl_spline *spline_x, gsl_interp_accel *spline_acc,
    real_t **x, real_t **xPrime, real_t **x_f, real_t **xPrime_f
) {
    *x        = new real_t[nr];
    *xPrime   = new real_t[nr];
    *x_f      = new real_t[nr+1];
    *xPrime_f = new real_t[nr+1];

    for (len_t ir=0; ir < nr; ir++){
        if (nProvided == 1) {
            (*x)[ir]      = xProvided[0];
            (*xPrime)[ir] = 0;
        } else {
            (*x)[ir]      = gsl_spline_eval(spline_x, r[ir], spline_acc);
            (*xPrime)[ir] = gsl_spline_eval_deriv(spline_x, r[ir], spline_acc);
        }
    }
    for (len_t ir=0; ir < nr+1; ir++){
        if (nProvided == 1) {
            (*x_f)[ir]      = xProvided[0];
            (*xPrime_f)[ir] = 0;
        } else {
            (*x_f)[ir]      = gsl_spline_eval(spline_x, r_f[ir], spline_acc);
            (*xPrime_f)[ir] = gsl_spline_eval_deriv(spline_x, r_f[ir], spline_acc);
        }
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


/**
 * Builds GSL spline objects from input shape profiles
 */
void AnalyticBRadialGridGenerator::constructSplines(struct shape_profiles *pp){
    auto construct_spline = [](const len_t n, const real_t *r, const real_t *x) {
        const gsl_interp_type *tp;
        if (n == 2)
            tp = gsl_interp_linear;
        else
            tp = gsl_interp_steffen;

        gsl_spline *s = gsl_spline_alloc(tp, n);
        gsl_spline_init(s, r, x, n);

        return s;
    };

    // Allocate splines for shape parameters (if necessary)
    if (pp->nG > 1) {
        this->spline_G = construct_spline(pp->nG, pp->G_r, pp->G);
        this->gsl_acc_G = gsl_interp_accel_alloc();
    }
    if (pp->npsi > 1) {
        this->spline_psi = construct_spline(pp->npsi, pp->psi_r, pp->psi);
        this->gsl_acc_psi = gsl_interp_accel_alloc();
    }
    if (pp->nkappa > 1) {
        this->spline_kappa = construct_spline(pp->nkappa, pp->kappa_r, pp->kappa);
        this->gsl_acc_kappa = gsl_interp_accel_alloc();
    }
    if (pp->ndelta > 1) {
        this->spline_delta = construct_spline(pp->ndelta, pp->delta_r, pp->delta);
        this->gsl_acc_delta = gsl_interp_accel_alloc();
    }
    if (pp->nDelta > 1) {
        this->spline_Delta = construct_spline(pp->nDelta, pp->Delta_r, pp->Delta);
        this->gsl_acc_Delta = gsl_interp_accel_alloc();
    }
}
