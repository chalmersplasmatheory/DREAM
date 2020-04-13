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
 *
 * nx: Number of radial grid points.
 * G: Toroidal magnetic field component as function of minor radius 
 * Psi_p0: Reference poloidal magnetic flux as function of minor radius 
 * x0: Value of inner radial flux grid point.
 * xa: Value of outer radial flux grid point.
 */

AnalyticBRadialGridGenerator::AnalyticBRadialGridGenerator(
     const len_t nr,  real_t r0,  real_t ra, real_t R0,
    real_t *rProfiles, len_t nrProfiles, real_t *Gs, real_t *psi_p0s,
             real_t *kappas, real_t *deltas, real_t *Deltas
) : RadialGridGenerator(nr), rMin(r0), rMax(ra), R0(R0){
    
    this->GsProvided     = Gs;
    this->psisProvided   = psi_p0s;
    this->kappasProvided = kappas;
    this->deltasProvided = deltas;
    this->DeltasProvided = Deltas;
    this->nrProfiles     = nrProfiles;
    this->rProfilesProvided = rProfiles;
    
    spline_x = gsl_spline_alloc(gsl_interp_steffen, nrProfiles);
    gsl_acc  = gsl_interp_accel_alloc(); 
}

AnalyticBRadialGridGenerator::~AnalyticBRadialGridGenerator(){
    gsl_spline_free (spline_x);
    gsl_interp_accel_free (gsl_acc);
}

/**
 * Interpolates input shape-parameter profiles (kappa, delta, ...) which are defined on 
 * input rProfilesProvided array to the r and r_f grids
 */
void AnalyticBRadialGridGenerator::InterpolateInputProfileToGrid(real_t *r, real_t *r_f, real_t *x,real_t *xPrime, real_t *x_f, real_t *xPrime_f,real_t *xProvided){
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
        xPrime_f[ir] = (gsl_spline_eval(spline_x,r[ir]+h,gsl_acc) - x_f[ir])/h;
        }
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

    InterpolateInputProfileToGrid(r,r_f,G,GPrime,G_f,GPrime_f,GsProvided);
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
 * Re-build the phase space jacobians.
 *
 * grid: Grid to build jacobians for.
 */
void AnalyticBRadialGridGenerator::RebuildJacobians(RadialGrid *rGrid, MomentumGrid **momentumGrids, MagneticQuantityHandler *mgnQtyHandler) {
    CreateMagneticFieldData(rGrid->GetR(),rGrid->GetR_f());

    rGrid->InitializeMagneticField(ntheta_ref, theta_ref,
            B_ref, B_ref_f,
            Bmin_ref, Bmin_ref_f,
            Bmax_ref, Bmax_ref_f
        );

    mgnQtyHandler->Initialize(momentumGrids,
                     ntheta_ref, theta_ref, 
                     B_ref, Jacobian_ref,
                     ROverR0_ref, NablaR2_ref,
                     B_ref_f, Jacobian_ref_f,
                     ROverR0_ref_f, NablaR2_ref_f);
    
    rGrid->InitializeVprime(mgnQtyHandler->GetVp(0),mgnQtyHandler->GetVp(1),
                            mgnQtyHandler->GetVp(2),mgnQtyHandler->GetVp(3),
                            mgnQtyHandler->GetVpVol(false), mgnQtyHandler->GetVpVol(true));
}


void AnalyticBRadialGridGenerator::CreateMagneticFieldData(const real_t *r, const real_t *r_f) {
    DeallocateMagneticFieldData();

    B_ref          = new real_t*[GetNr()];
    Jacobian_ref   = new real_t*[GetNr()];
    ROverR0_ref    = new real_t*[GetNr()];
    NablaR2_ref    = new real_t*[GetNr()];
    Bmin_ref       = new real_t[GetNr()];
    Bmax_ref       = new real_t[GetNr()];
    B_ref_f        = new real_t*[(GetNr()+1)];
    Jacobian_ref_f = new real_t*[(GetNr()+1)];
    ROverR0_ref_f  = new real_t*[(GetNr()+1)];
    NablaR2_ref_f  = new real_t*[(GetNr()+1)];
    Bmin_ref_f     = new real_t[GetNr()+1];
    Bmax_ref_f     = new real_t[GetNr()+1];
    

    theta_ref = new real_t[ntheta_ref];
    real_t dth = 2*M_PI / (ntheta_ref-1);
    for(len_t it=0; it<ntheta_ref; it++) {
        theta_ref[it] = it*dth; 
    }
    real_t R, st, ct;
    for (len_t ir = 0; ir < GetNr(); ir++){
        Jacobian_ref[ir] = new real_t[ntheta_ref];
        B_ref[ir]        = new real_t[ntheta_ref];
        ROverR0_ref[ir]  = new real_t[ntheta_ref];
        NablaR2_ref[ir]  = new real_t[ntheta_ref];
        for(len_t it=0; it<ntheta_ref; it++){
            ct = cos(theta_ref[it]);
            st = sin(theta_ref[it]);;
            R = R0 + Delta[ir] + r[ir]*cos(theta_ref[it] + delta[ir]*st);
            ROverR0_ref[ir][it] = R;
            Jacobian_ref[ir][it] = kappa[ir]*r[ir]*R * ( cos(delta[ir]*st) + DeltaPrime[ir]*ct
            + st*sin(theta_ref[it]+delta[ir]*st) * ( r[ir]*kappaPrime[ir]/kappa[ir] + delta[ir]*ct
            * ( 1 + r[ir]*kappaPrime[ir]/kappa[ir] - r[ir]*deltaPrime[ir]/delta[ir] ) ) );
            
            NablaR2_ref[ir][it] = kappa[ir]*kappa[ir]*r[ir]*r[ir]*R*R
                *( ct*ct + (1+delta[ir]*delta[ir])*(1+delta[ir]*delta[ir])/(kappa[ir]*kappa[ir]) 
                * sin(theta_ref[it]+delta[ir]*st)*sin(theta_ref[it]+delta[ir]*st) ) 
                / ( Jacobian_ref[ir][it]*Jacobian_ref[ir][it]); 
            
            B_ref[ir][it] = G[ir]*G[ir]/(R*R)
                                + NablaR2_ref[ir][it] * psiPrime[ir]*psiPrime[ir];
        }
        Bmin_ref[ir] = B_ref[ir][0];
        Bmax_ref[ir] = B_ref[ir][0];
        for(len_t it=0; it<ntheta_ref; it++){
            if (Bmin_ref[ir] > B_ref[ir][it])
                Bmin_ref[ir] = B_ref[ir][it];
            if (Bmax_ref[ir] <= B_ref[ir][it])
                Bmax_ref[ir] = B_ref[ir][it];
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
            R = R0 + Delta_f[ir] + r_f[ir]*cos(theta_ref[it] + delta_f[ir]*st);
            ROverR0_ref_f[ir][it] = R;
            Jacobian_ref_f[ir][it] = kappa_f[ir]*r_f[ir]*R * ( cos(delta_f[ir]*st) + DeltaPrime_f[ir]*ct
            + st*sin(theta_ref[it]+delta_f[ir]*st) * ( r_f[ir]*kappaPrime_f[ir]/kappa_f[ir] + delta_f[ir]*ct
            * ( 1 + r_f[ir]*kappaPrime_f[ir]/kappa_f[ir] - r_f[ir]*deltaPrime_f[ir]/delta_f[ir] ) ) );
            
            NablaR2_ref_f[ir][it] = kappa_f[ir]*kappa_f[ir]*r_f[ir]*r_f[ir]*R*R
                *( ct*ct + (1+delta_f[ir]*delta_f[ir])*(1+delta_f[ir]*delta_f[ir])/(kappa_f[ir]*kappa_f[ir]) 
                * sin(theta_ref[it]+delta_f[ir]*st)*sin(theta_ref[it]+delta_f[ir]*st) ) 
                / ( Jacobian_ref_f[ir][it]*Jacobian_ref_f[ir][it]); 
            
            B_ref_f[ir][it] = G_f[ir]*G_f[ir]/(R*R)
                                + NablaR2_ref_f[ir][it] * psiPrime_f[ir]*psiPrime_f[ir];
        }    
        Bmin_ref_f[ir] = B_ref_f[ir][0];
        Bmax_ref_f[ir] = B_ref_f[ir][0];
        for(len_t it=0; it<ntheta_ref; it++){
            if (Bmin_ref_f[ir] > B_ref_f[ir][it])
                Bmin_ref_f[ir] = B_ref_f[ir][it];
            if (Bmax_ref_f[ir] <= B_ref_f[ir][it])
                Bmax_ref_f[ir] = B_ref_f[ir][it];
        }
    }

}



void AnalyticBRadialGridGenerator::DeallocateMagneticFieldData(){
    if (B_ref==nullptr)
        return;

    for(len_t ir = 0; ir<GetNr(); ir++){
        delete [] B_ref[ir];
        delete [] Jacobian_ref[ir];
        delete [] ROverR0_ref[ir];
        delete [] NablaR2_ref[ir];
    }
    for(len_t ir = 0; ir<GetNr()+1; ir++){
        delete [] B_ref_f[ir];
        delete [] Jacobian_ref_f[ir];
        delete [] ROverR0_ref_f[ir];
        delete [] NablaR2_ref_f[ir];
    }
    delete [] theta_ref;
    delete [] Bmin_ref;
    delete [] Bmin_ref_f;
    delete [] Bmax_ref;
    delete [] Bmax_ref_f;
    delete [] B_ref;
    delete [] Jacobian_ref;
    delete [] ROverR0_ref;
    delete [] NablaR2_ref;
    delete [] B_ref_f;
    delete [] Jacobian_ref_f;
    delete [] ROverR0_ref_f;
    delete [] NablaR2_ref_f;
}

