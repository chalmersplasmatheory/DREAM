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


using namespace DREAM::FVM;
using namespace std;

/**
 * Constructor.
 * nr: Number of radial grid points.
 * GOverR0: Toroidal magnetic field component as function of minor radius 
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

    this->R0IsInf = isinf(R0);
    constructSplines(profiles);
    isUpDownSymmetric = true;
}

/**
 * Constructor.
 * 
 * r_f_input: Grid points on the radial flux grid (e.g. the cell edges)
 * nr: Number of radial grid points.
 * GOverR0: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 */
AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const real_t *r_f_input, const len_t nr, real_t R0, len_t ntheta_interp,
     struct shape_profiles *profiles
) : RadialGridGenerator(nr), rMin(r_f_input[0]), rMax(r_f_input[nr]), providedProfiles(profiles) {
    this->R0             = R0;
    this->ntheta_interp  = ntheta_interp;

    this->R0IsInf = isinf(R0);
    this->rf_provided = new real_t[nr+1];
    for(len_t i=0; i<nr+1; i++)
        this->rf_provided[i] = r_f_input[i];

    //delete [] r_f_input;

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
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nG,     pp->G_r,     pp->GOverR0, spline_G,     gsl_acc_G,     &BtorGOverR0, &GPrime,      &BtorGOverR0_f, &GPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->npsi,   pp->psi_r,   pp->psi,     spline_psi,   gsl_acc_psi,   &psi,         &psiPrimeRef, &psi_f,         &psiPrimeRef_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nkappa, pp->kappa_r, pp->kappa,   spline_kappa, gsl_acc_kappa, &kappa,       &kappaPrime,  &kappa_f,       &kappaPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->ndelta, pp->delta_r, pp->delta,   spline_delta, gsl_acc_delta, &delta,       &deltaPrime,  &delta_f,       &deltaPrime_f);
    InterpolateInputProfileToGrid(GetNr(), r, r_f, pp->nDelta, pp->Delta_r, pp->Delta,   spline_Delta, gsl_acc_Delta, &Delta,       &DeltaPrime,  &Delta_f,       &DeltaPrime_f);

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
    real_t ct = cos(theta);
    real_t st = sin(theta);

    return r[ir]*ROverR0AtTheta(ir,theta) * normalizedJacobian(ir,theta,ct,st);
}
// Same as JacobianAtTheta but evaluated on the radial flux grid
real_t AnalyticBRadialGridGenerator::JacobianAtTheta_f(const len_t ir, const real_t theta){
    real_t ct = cos(theta);
    real_t st = sin(theta);

    return r_f[ir]*ROverR0AtTheta_f(ir,theta) * normalizedJacobian_f(ir,theta,ct,st);
}

/**
 * Evaluates |nabla r|^2 at radial grid point ir and poloidal angle theta
 */
real_t AnalyticBRadialGridGenerator::NablaR2AtTheta(const len_t ir, const real_t theta){
    real_t ct = cos(theta);
    real_t st = sin(theta);

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
    real_t ct = cos(theta);
    real_t st = sin(theta);

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
    real_t ct = cos(theta);
    real_t st = sin(theta);

    real_t sdt = 0.0;
    real_t cdt = 1.0;
    if(delta[ir]) {
        cdt = cos(delta[ir]*st);
        sdt = sin(delta[ir]*st);
    }

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
    real_t ct = cos(theta);
    real_t st = sin(theta);

    real_t sdt = 0.0;
    real_t cdt = 1.0;
    if(delta_f[ir]) {
        cdt = cos(delta_f[ir]*st);
        sdt = sin(delta_f[ir]*st);
    }

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
    const len_t nProvided, const real_t *rProvided, const real_t *xProvided,
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
            if (r[ir] > rProvided[nProvided-1]) {
                (*x)[ir]      = xProvided[nProvided-1];
                (*xPrime)[ir] = gsl_spline_eval_deriv(spline_x, rProvided[nProvided-1], spline_acc);
			} else if (r[ir] < rProvided[0]) {
				// Extrapolate linearly
				real_t dx     = gsl_spline_eval_deriv(spline_x, rProvided[0], spline_acc);
				(*x)[ir]      = xProvided[0] - (rProvided[0]-r[ir])*dx;
				(*xPrime)[ir] = dx;
            } else {
                (*x)[ir]      = gsl_spline_eval(spline_x, r[ir], spline_acc);
                (*xPrime)[ir] = gsl_spline_eval_deriv(spline_x, r[ir], spline_acc);
            }
        }
    }
    for (len_t ir=0; ir < nr+1; ir++){
        if (nProvided == 1) {
            (*x_f)[ir]      = xProvided[0];
            (*xPrime_f)[ir] = 0;
        } else {
            if (r_f[ir] > rProvided[nProvided-1]) {
                (*x_f)[ir]      = xProvided[nProvided-1];
                (*xPrime_f)[ir] = gsl_spline_eval_deriv(spline_x, rProvided[nProvided-1], spline_acc);
			} else if (r_f[ir] < rProvided[0]) {
				real_t dx       = gsl_spline_eval_deriv(spline_x, rProvided[0], spline_acc);
				(*x_f)[ir]      = xProvided[0] - (rProvided[0]-r_f[ir])*dx;
				(*xPrime)[ir]   = dx;
            } else {
                (*x_f)[ir]      = gsl_spline_eval(spline_x, r_f[ir], spline_acc);
                (*xPrime_f)[ir] = gsl_spline_eval_deriv(spline_x, r_f[ir], spline_acc);
            }
        }
    }
}

/**
 * Interpolates the input chaping profile at a single radius r, and extrapolates using 
 * the maximum and minimum r for which the profile is defined. This is useful for the 
 * GetRThetaPhiFromCartesian() and GetGradRCartesian() functions used for e.g. the SPI
 * module in DREAM, as these functions must be able to handle coordinates outside the plasma
 *
 */
real_t AnalyticBRadialGridGenerator::InterpolateInputProfileSingleExtrap(real_t r,
    const len_t nProvided, const real_t *xProvided, const real_t *xProvided_r,
    gsl_spline *spline_x, gsl_interp_accel *spline_acc
) {
    if(nProvided==1){
        return xProvided[0];
    } else {
        if(r>xProvided_r[nProvided-1])
            return xProvided_r[nProvided-1];
        else if(r<xProvided_r[0])
            return xProvided_r[0];
        else
            return gsl_spline_eval(spline_x, r, spline_acc);
    }   
}

/**
 * Similar to InterpolateInputProfileSingleExtrap(), but for derivatives
 */
real_t AnalyticBRadialGridGenerator::InterpolateInputProfileSingleDerivExtrap(real_t r,
    const len_t nProvided, const real_t *xProvided_r,
    gsl_spline *spline_x, gsl_interp_accel *spline_acc
) {
    if(nProvided==1){
        return 0;
    } else {
        if(r>xProvided_r[nProvided-1])
            return gsl_spline_eval_deriv(spline_x, xProvided_r[nProvided-1], spline_acc);
        else if(r<xProvided_r[0])
            return gsl_spline_eval_deriv(spline_x, xProvided_r[0], spline_acc);
        else
            return gsl_spline_eval_deriv(spline_x, r, spline_acc);
    }   
}

/**
* Helper function for interpolation of elongation data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputElongation(real_t r){
	const real_t *kappaProf = this->providedProfiles->kappa;
	const real_t *kappa_r = this->providedProfiles->kappa_r;
	const len_t nkappa = this->providedProfiles->nkappa;
	
	return InterpolateInputProfileSingleExtrap(r, nkappa, kappaProf, kappa_r, 
	            this->spline_kappa, gsl_acc_kappa);
}

/**
* Helper function for interpolation of the derivative of the elongation data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputElongationDeriv(real_t r){
	const real_t *kappa_r = this->providedProfiles->kappa_r;
	const len_t nkappa = this->providedProfiles->nkappa;
	
	return InterpolateInputProfileSingleDerivExtrap(r, nkappa, kappa_r, 
	            this->spline_kappa, gsl_acc_kappa);
}

/**
* Helper function for interpolation of triangularity data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputTriangularity(real_t r){
	const real_t *deltaProf = this->providedProfiles->delta;
	const real_t *delta_r = this->providedProfiles->delta_r;
	const len_t ndelta = this->providedProfiles->ndelta;
	
	return InterpolateInputProfileSingleExtrap(r, ndelta, deltaProf, delta_r, 
	            this->spline_delta, gsl_acc_delta);
}

/**
* Helper function for interpolation of the derivative of the triangularity data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputTriangularityDeriv(real_t r){
	const real_t *delta_r = this->providedProfiles->delta_r;
	const len_t ndelta = this->providedProfiles->ndelta;
	
	return InterpolateInputProfileSingleDerivExtrap(r, ndelta, delta_r, 
	            this->spline_delta, gsl_acc_delta);
}

/**
* Helper function for interpolation of Shafranov shift data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputShafranovShift(real_t r){
	const real_t *DeltaProf = this->providedProfiles->Delta;
	const real_t *Delta_r = this->providedProfiles->Delta_r;
	const len_t nDelta = this->providedProfiles->nDelta;
	
	return InterpolateInputProfileSingleExtrap(r, nDelta, DeltaProf, Delta_r, 
	            this->spline_Delta, gsl_acc_Delta);
}

/**
* Helper function for interpolation of the derivative of the Shafranov shift data
*/
real_t AnalyticBRadialGridGenerator::InterpolateInputShafranovShiftDeriv(real_t r){
	const real_t *Delta_r = this->providedProfiles->Delta_r;
	const len_t nDelta = this->providedProfiles->nDelta;
	
	return InterpolateInputProfileSingleDerivExtrap(r, nDelta, Delta_r, 
	            this->spline_Delta, gsl_acc_Delta);
}

/**
 * Calculate the flux surface coordinates 'r, theta and phi' corresponding to the given
 * Cartesian coordinates (x,y,z).
 *
 * (The Cartesian coordinate system is oriented such that x and y span
 * the poloidal plane. The origin of x and y is the magnetic axis.)
 *
 * The coordinate transformation is solved numerically using a bisection scheme until
 * a precision of lengthScale*CartesianCoordinatesTol is achieved. The initial search 
 * intervall is r = startingGuess +/- lengthScale, translated by a multiple of 2*lengthScale
 * until the search interval contains a sign change
 *
 * See doc/notes/SPICoordinates.pdf for more information
 */
void AnalyticBRadialGridGenerator::GetRThetaPhiFromCartesian(real_t* r, real_t* theta, real_t* phi,
    real_t x, real_t y, real_t z, real_t lengthScale, real_t startingGuessR
) {

    real_t RMinusR0;
    if(R0IsInf)
        RMinusR0 = x;
    else
        RMinusR0 = hypot(R0+x, z)-R0;
	
                
    // We start by finding the flux surface label 'rmin' corresponding to cos(theta)=0 (using Newton iteration)
    // and use it to determine in which quadrant we should look for the solution
	real_t r_tmp=std::abs(y);
	real_t Zdiff;// y in SPI coordinates is Z in cylindrical coordinates
	real_t kappa;
	real_t kappa_p;
	
	do{
	    kappa = InterpolateInputElongation(r_tmp);
	    kappa_p = InterpolateInputElongationDeriv(r_tmp);
	    Zdiff = std::abs(y)-r_tmp*kappa;
	    r_tmp = r_tmp + Zdiff/(r_tmp*kappa_p+kappa);
	    
	    	
	}while (std::abs(Zdiff)>lengthScale*CartesianCoordinateTol);
	
	real_t rmin=r_tmp;
	
	// Evaluate RMinusR0_crit corresponding to rmin
	real_t delta;
	real_t Delta; 
    delta = InterpolateInputTriangularity(rmin);
    Delta = InterpolateInputShafranovShift(rmin);
    
    real_t RMinusR0_crit=Delta-rmin*sin(delta);
	
	// Determine the quadrant
	len_t quadrant=0;
	if(RMinusR0>RMinusR0_crit && y>=0)
	    quadrant=1;
	else if(RMinusR0<RMinusR0_crit && y>=0)
	    quadrant=2;
    else if(RMinusR0<RMinusR0_crit && y<0)
        quadrant=3;
    else if(RMinusR0>RMinusR0_crit && y<0)
        quadrant=4;
	
	// Bisection to find radial coordinate corresponding
	// to 'r' at 'Z=y'...
	// We make a guess for a valid search intervall of startingGuessR+/-lengthScale, 
	// and check if it has to be moved before actually starting with the bisection
	real_t ra = std::max(startingGuessR-lengthScale,rmin), rb=ra+2*lengthScale;
	real_t RMinusR0a, RMinusR0b;
	real_t st;
	real_t ct;
	do{
	    // Calculate RMinusR0a corresponding to ra
        delta = InterpolateInputTriangularity(ra);
        Delta = InterpolateInputShafranovShift(ra);
	    kappa = InterpolateInputElongation(ra);
        
        if(ra==0){
        	st=1; // The value does not matter here as RMinusR0a = Delta anyway
        }else{
        	st = y/(ra*kappa);
        	
        	// As we have already made sure that ra>=rmin, an st outside [-1,1] 
        	// should only depend on round-off errors, so it should be safe to do this
        	if(st>1.0)
        		st=1.0;
        	else if(st<-1.0)
        		st=-1.0;
        }
        	
	    if(quadrant==1 || quadrant==4)
	        ct = sqrt(1-st*st);
	    else
	        ct = -sqrt(1-st*st);
	        
        RMinusR0a = Delta + ra*(ct*cos(delta*st)-st*sin(delta*st));
        
        // Similarly for rb
        delta = InterpolateInputTriangularity(rb);
        Delta = InterpolateInputShafranovShift(rb);
	    kappa = InterpolateInputElongation(rb);
        
        st = y/(rb*kappa);
	    if(quadrant==1 || quadrant==4)
	        ct = sqrt(1-st*st);
	    else
	        ct = -sqrt(1-st*st);
	        
        RMinusR0b = Delta + rb*(ct*cos(delta*st)-st*sin(delta*st));
        
        // Move the search intervall if necessary (but make sure to keep ra>=rmin)
        if(quadrant==1 || quadrant==4){
            if(RMinusR0a>RMinusR0 && RMinusR0b>RMinusR0){
                ra=std::max(rmin, ra-2*lengthScale);
                rb=ra+2*lengthScale;
		    }
	        else if(RMinusR0a<RMinusR0 && RMinusR0b<RMinusR0){
	            rb+=2*lengthScale;
	            ra+=2*lengthScale;
            }
        }else{
            if(RMinusR0a>RMinusR0 && RMinusR0b>RMinusR0){
	            rb+=2*lengthScale;
	            ra+=2*lengthScale;
		    }
	        else if(RMinusR0a<RMinusR0 && RMinusR0b<RMinusR0){
                ra-=2*lengthScale;
                rb-=2*lengthScale;
            }
        }
        
	}while((RMinusR0a>RMinusR0 && RMinusR0b>RMinusR0) || (RMinusR0a<RMinusR0 && RMinusR0b<RMinusR0));
	
	// Make the bisection
	real_t RMinusR0_tmp;
	do{
	    r_tmp=(ra+rb)/2.0;
        delta = InterpolateInputTriangularity(r_tmp);
        Delta = InterpolateInputShafranovShift(r_tmp);
	    kappa = InterpolateInputElongation(r_tmp);
        
        st = y/(r_tmp*kappa);
	    if(quadrant==1 || quadrant==4)
	        ct = sqrt(1-st*st);
	    else
	        ct = -sqrt(1-st*st);
	        
	    RMinusR0_tmp = Delta + r_tmp*(ct*cos(delta*st)-st*sin(delta*st));
        if(quadrant==1 || quadrant==4){
            if(RMinusR0_tmp>RMinusR0){
                rb=r_tmp;
		    }
	        else{
	            ra=r_tmp;
            }
        }else{
            if(RMinusR0_tmp>RMinusR0){
                ra=r_tmp;
		    }
	        else{
                rb=r_tmp;
            }
        }
	}while(std::abs(rb-ra) > lengthScale*CartesianCoordinateTol);
	
	// Newton solver (disabled for now)
	/*r_tmp=startingGuessR;
	
	real_t Delta_p;
	real_t delta_p;
	real_t st;
	real_t st_p;
	real_t ct;
	real_t ct_p;
	real_t RMinusR0_newton;
	real_t ddrRMinusR0_newton;
	do {
	    kappa   = InterpolateInputElongation(r_tmp);
	    kappa_p = InterpolateInputElongationDeriv(r_tmp);
	    delta   = InterpolateInputTriangularity(r_tmp);
	    delta_p = InterpolateInputTriangularityDeriv(r_tmp);
	    Delta   = InterpolateInputShafranovShift(r_tmp);
	    Delta_p = InterpolateInputShafranovShiftDeriv(r_tmp);
	    
	    st = y/(r_tmp*kappa);
	    st_p = -y/(r_tmp*r_tmp*kappa*kappa)*(kappa+r*kappa_p);
	    
	    if(quadrant==1 || quadrant==4)
	        ct = sqrt(1-st*st);
	    else
	        ct = -sqrt(1-st*st);
	    ct_p = -st/ct*st_p;
	    
        RMinusR0_newton = Delta + r_tmp*(ct*cos(delta*st)-st*sin(delta*st));
        ddrRMinusR0_newton = Delta_p + ct*cos(delta*st)-st*sin(delta*st)+
                             r_tmp*(ct_p*cos(delta*st) - st_p*sin(delta*st)-
                             (ct*sin(delta*st)+st*cos(delta*st))*(delta_p*st+delta*st_p));
        r_tmp = r_tmp - (RMinusR0_newton - RMinusR0)/ddrRMinusR0_newton;
	} while(std::abs(RMinusR0_newton-RMinusR0) > lengthScale * CartesianCoordinateTol);*/
	
	// Set the output values for the flux surface coordinates
	*r=r_tmp;
	
	// We define theta to be in the range 0<=theta<2*pi
	if(quadrant==1)
	    *theta=asin(st);
	else if(quadrant==4)
	    *theta=2*M_PI+asin(st);
	else if(RMinusR0<Delta-r_tmp*sin(delta))
	    *theta=M_PI-asin(st);
	    
	*phi = atan2(z,(R0+x)); 
    
}

/**
 * Calculates the gradient of the minor radius coordinate 'r' in cartesian coordinates
 *
 * See section 1.2 in doc/notes/theory.pdf for derivation
 */
void AnalyticBRadialGridGenerator::GetGradRCartesian(real_t* gradr, real_t r, real_t theta, real_t phi) {
	
    real_t delta  = InterpolateInputTriangularity(r);
    real_t deltap = InterpolateInputTriangularityDeriv(r);
    
    real_t Deltap = InterpolateInputShafranovShiftDeriv(r);
    
    real_t kappa  = InterpolateInputElongation(r);
    real_t kappap = InterpolateInputElongationDeriv(r);
    
    real_t st = sin(theta);
    real_t ct = cos(theta);
    
    real_t dRdr = Deltap + cos(theta + delta*st) - r*deltap*st*sin(theta + delta*st);
    real_t dRdtheta = -r*(1+delta*ct)*sin(theta + delta*st);
    real_t dzdr = kappa*(1+r*kappap/kappa)*st;
    real_t dzdtheta = r*kappa*ct;
    
    real_t prefactor = 1/(dRdr*dzdtheta - dRdtheta*dzdr);
    real_t Rhatfactor = prefactor*dzdtheta;
    
    gradr[0] = cos(phi)*Rhatfactor;
    gradr[1] = -dRdtheta*prefactor;
    gradr[2] = sin(phi)*Rhatfactor; 
    
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
        this->spline_G = construct_spline(pp->nG, pp->G_r, pp->GOverR0);
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
