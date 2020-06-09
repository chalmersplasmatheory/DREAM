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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

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
     const len_t nr,  real_t r0,  real_t ra, real_t R0, len_t ntheta_ref, len_t ntheta_interp,
    real_t *rProfiles, len_t nrProfiles, real_t *Gs, real_t *psi_p0s,
             real_t *kappas, real_t *deltas, real_t *Deltas
) : RadialGridGenerator(nr), rMin(r0), rMax(ra), R0(R0) {
    this->ntheta_ref     = ntheta_ref;
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
    DeallocateMagneticFieldData();
    DeallocateShapeProfiles();
}

/**
 * (Re-)builds the given radial grid.
 *
 * rGrid: Radial grid to re-build.
 */
bool AnalyticBRadialGridGenerator::Rebuild(const real_t, RadialGrid *rGrid) {
    
    real_t
        *r    = new real_t[GetNr()],
        *r_f  = new real_t[GetNr()+1],
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


    this->isBuilt = true;

    DeallocateShapeProfiles();
    InterpolateInputProfileToGrid(r,r_f,BtorGOverR0,GPrime,BtorGOverR0_f,GPrime_f,GsProvided);
    InterpolateInputProfileToGrid(r,r_f,psi,psiPrime,psi_f,psiPrime_f,psisProvided);
    InterpolateInputProfileToGrid(r,r_f,kappa,kappaPrime,kappa_f,kappaPrime_f,kappasProvided);
    InterpolateInputProfileToGrid(r,r_f,delta,deltaPrime,delta_f,deltaPrime_f,deltasProvided);
    InterpolateInputProfileToGrid(r,r_f,Delta,DeltaPrime,Delta_f,DeltaPrime_f,DeltasProvided);
    
    return true;
}

real_t AnalyticBRadialGridGenerator::diffFunc(real_t r, std::function<real_t(real_t)> F){
    real_t sqrteps = sqrt(__DBL_EPSILON__);
    real_t h = sqrteps * ( 1 + abs(r) ); 
    return (F(r+h/2)-F(r-h/2))/h;
}


/**
 * The method generates data on a reference theta grid which is 
 * uniformly spaced from 0 to 2pi with ntheta_ref points. Data 
 * is created for quantities which appear in bounce averages etc:
 * the magnetic field B, spatial Jacobian J (normalized to R0), the gradient
 * |nabla r|^2 NablaR2 and normalized major radius R/R0, ROverR0.
 */
void AnalyticBRadialGridGenerator::CreateMagneticFieldData(const real_t *r, const real_t *r_f) {

    B_ref          = new real_t*[GetNr()];
    Jacobian_ref   = new real_t*[GetNr()];
    ROverR0_ref    = new real_t*[GetNr()];
    NablaR2_ref    = new real_t*[GetNr()];
    Bmin           = new real_t[GetNr()];
    Bmax           = new real_t[GetNr()];
//    BtorGOverR0    = new real_t[GetNr()];
    B_ref_f        = new real_t*[(GetNr()+1)];
    Jacobian_ref_f = new real_t*[(GetNr()+1)];
    ROverR0_ref_f  = new real_t*[(GetNr()+1)];
    NablaR2_ref_f  = new real_t*[(GetNr()+1)];
    Bmin_f         = new real_t[GetNr()+1];
    Bmax_f         = new real_t[GetNr()+1];
//    BtorGOverR0_f  = new real_t[GetNr()+1];

    theta_ref = new real_t[ntheta_ref];
    real_t dth = 2*M_PI / (ntheta_ref-1);
    for(len_t it=0; it<ntheta_ref; it++) {
        theta_ref[it] = it*dth; 
    }
    real_t st, ct, JOverRr;
    for (len_t ir = 0; ir < GetNr(); ir++){
        Jacobian_ref[ir] = new real_t[ntheta_ref];
        B_ref[ir]        = new real_t[ntheta_ref];
        ROverR0_ref[ir]  = new real_t[ntheta_ref];
        NablaR2_ref[ir]  = new real_t[ntheta_ref];
        for(len_t it=0; it<ntheta_ref; it++){
            ct = cos(theta_ref[it]);
            st = sin(theta_ref[it]);;
            if(isinf(R0))
                ROverR0_ref[ir][it] = 1; 
            else 
                ROverR0_ref[ir][it] = 1 + (Delta[ir] + r[ir]*cos(theta_ref[it] + delta[ir]*st))/R0;
            JOverRr = ( kappa[ir]*cos(delta[ir]*st) + kappa[ir]*DeltaPrime[ir]*ct
            + st*sin(theta_ref[it]+delta[ir]*st) * ( r[ir]*kappaPrime[ir] + 
            ct * ( delta[ir]*kappa[ir] + r[ir]*delta[ir]*kappaPrime[ir] - kappa[ir]*r[ir]*deltaPrime[ir] ) ) );
            Jacobian_ref[ir][it] = r[ir]*ROverR0_ref[ir][it] * JOverRr;
            
            NablaR2_ref[ir][it] = ( kappa[ir]*kappa[ir]* ct*ct + (1+delta[ir]*ct)*(1+delta[ir]*ct) 
                * sin(theta_ref[it]+delta[ir]*st)*sin(theta_ref[it]+delta[ir]*st) ) 
                / ( JOverRr*JOverRr); 
            
            B_ref[ir][it] = sqrt(BtorGOverR0[ir]*BtorGOverR0[ir] + NablaR2_ref[ir][it] * psiPrime[ir]*psiPrime[ir]/(R0*R0)) / ROverR0_ref[ir][it];
        }
        Bmin[ir] = B_ref[ir][0];
        Bmax[ir] = B_ref[ir][0];
//        BtorGOverR0[ir] = G[ir];

        for(len_t it=0; it<ntheta_ref; it++){
            if (Bmin[ir] > B_ref[ir][it])
                Bmin[ir] = B_ref[ir][it];
            if (Bmax[ir] <= B_ref[ir][it])
                Bmax[ir] = B_ref[ir][it];
        }
    }

    for (len_t ir = 0; ir < GetNr()+1; ir++){
        Jacobian_ref_f[ir] = new real_t[ntheta_ref];
        B_ref_f[ir]        = new real_t[ntheta_ref];
        ROverR0_ref_f[ir]  = new real_t[ntheta_ref];
        NablaR2_ref_f[ir]  = new real_t[ntheta_ref];
        
        for(len_t it=0; it<ntheta_ref; it++){
            ct = cos(theta_ref[it]);
            st = sin(theta_ref[it]);;
            if(isinf(R0))
                ROverR0_ref_f[ir][it] = 1; 
            else 
                ROverR0_ref_f[ir][it] = 1 + (Delta_f[ir] + r_f[ir]*cos(theta_ref[it] + delta_f[ir]*st))/R0;
            
            JOverRr =  ( kappa_f[ir]*cos(delta_f[ir]*st) + kappa_f[ir]*DeltaPrime_f[ir]*ct
            + st*sin(theta_ref[it]+delta_f[ir]*st) * ( r_f[ir]*kappaPrime_f[ir] +
            ct * (  delta_f[ir]*kappa_f[ir] + r_f[ir]* delta_f[ir]*kappaPrime_f[ir] - r_f[ir]*kappa_f[ir]*deltaPrime_f[ir] ) ) );
            Jacobian_ref_f[ir][it] = r_f[ir] * ROverR0_ref_f[ir][it] * JOverRr;
            
            NablaR2_ref_f[ir][it] = 
                (kappa_f[ir]*kappa_f[ir] * ct*ct + (1+delta_f[ir]*ct)*(1+delta_f[ir]*ct) 
                * sin(theta_ref[it]+delta_f[ir]*st)*sin(theta_ref[it]+delta_f[ir]*st) )  / (JOverRr*JOverRr); 
            
            B_ref_f[ir][it] = sqrt(BtorGOverR0_f[ir]*BtorGOverR0_f[ir] + NablaR2_ref_f[ir][it] * psiPrime_f[ir]*psiPrime_f[ir]/(R0*R0)) / ROverR0_ref_f[ir][it];
        }    
        Bmin_f[ir] = B_ref_f[ir][0];
        Bmax_f[ir] = B_ref_f[ir][0];

        for(len_t it=0; it<ntheta_ref; it++){
            if (Bmin_f[ir]  > B_ref_f[ir][it])
                Bmin_f[ir]  = B_ref_f[ir][it];
            if (Bmax_f[ir] <= B_ref_f[ir][it]){
                Bmax_f[ir]  = B_ref_f[ir][it];

            }
        }
    }

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


void AnalyticBRadialGridGenerator::DeallocateShapeProfiles(){
    if (psi==nullptr)
        return;

//    delete [] G;
    delete [] psi;
    delete [] kappa;
    delete [] delta;
    delete [] Delta;
//    delete [] G_f;
    delete [] psi_f;
    delete [] kappa_f;
    delete [] delta_f;
    delete [] Delta_f;
    delete [] GPrime;
    delete [] psiPrime;
    delete [] kappaPrime;
    delete [] deltaPrime;
    delete [] DeltaPrime;
    delete [] GPrime_f;
    delete [] psiPrime_f;
    delete [] kappaPrime_f;
    delete [] deltaPrime_f;
    delete [] DeltaPrime_f;
    

}
