/**
 * Implementation of the FluxSurfaceAverager class which 
 * handles everything needed to carry out flux surface 
 * averages in DREAM. Is initialized by a RadialGridGenerator
 * via the method SetReferenceMagneticFieldData.
 */

#include "FVM/Grid/FluxSurfaceQuantity.hpp"

using namespace std;
using namespace DREAM::FVM;

FluxSurfaceQuantity::FluxSurfaceQuantity(
    RadialGrid *rGrid,
    const gsl_interp_type *interpType
) : rGrid(rGrid), interpolationMethod(interpType) { 
    gsl_acc = gsl_interp_accel_alloc();
}

FluxSurfaceQuantity::~FluxSurfaceQuantity(){
    gsl_interp_accel_free(gsl_acc);
    DeallocateReferenceData();
    DeallocateInterpolatedData();
}


void FluxSurfaceQuantity::InterpolateMagneticDataToTheta(real_t *theta, len_t ntheta_interp){
    DeallocateInterpolatedData();
    quantityData = new real_t*[nr];
    for(len_t ir = 0; ir<nr; ir++){
        quantityData[ir] = new real_t[ntheta_interp];
        for(len_t it=0; it<ntheta_interp; it++){
            quantityData[ir][it] = gsl_spline_eval(quantitySpline[ir], theta[it], gsl_acc);
        }
    }
    quantityData_fr = new real_t*[nr+1];
    for(len_t ir = 0; ir<=nr; ir++){
        quantityData_fr[ir] = new real_t[ntheta_interp];
        for(len_t it=0; it<ntheta_interp; it++){
            quantityData_fr[ir][it] = gsl_spline_eval(quantitySpline_fr[ir], theta[it], gsl_acc);
        }
    }
}

void FluxSurfaceQuantity::DeallocateInterpolatedData(){
    if(quantityData == nullptr)
        return;

    for(len_t ir = 0; ir<nr; ir++)
        delete [] quantityData[ir];
    for(len_t ir = 0; ir<=nr; ir++)
        delete [] quantityData_fr[ir];
    delete [] quantityData;
    delete [] quantityData_fr;    
}


void FluxSurfaceQuantity::Initialize(real_t **referenceData, real_t **referenceData_fr, real_t *theta_ref, len_t ntheta_ref){
    DeallocateReferenceData();
    this->nr = rGrid->GetNr();
    this->theta_ref  = theta_ref;
    this->ntheta_ref = ntheta_ref;
    
    theta_ref_min = theta_ref[0];
    theta_ref_max = theta_ref[0];
    for(len_t it=0; it<ntheta_ref; it++){
        if(theta_ref[it]>theta_ref_max)
            theta_ref_max = theta_ref[it];
        if(theta_ref[it]<theta_ref_min)
            theta_ref_min = theta_ref[it];
    }

    this->referenceData    = referenceData;
    this->referenceData_fr = referenceData_fr;

    quantitySpline    = new gsl_spline*[nr];
    quantitySpline_fr = new gsl_spline*[nr+1];
    
    for ( len_t ir = 0; ir<nr; ir++){
        quantitySpline[ir] = gsl_spline_alloc(interpolationMethod, ntheta_ref);
        gsl_spline_init(quantitySpline[ir], theta_ref, referenceData[ir], ntheta_ref);
    }
    for ( len_t ir = 0; ir<=nr; ir++){
        quantitySpline_fr[ir] = gsl_spline_alloc(interpolationMethod, ntheta_ref);
        gsl_spline_init(quantitySpline_fr[ir], theta_ref, referenceData_fr[ir], ntheta_ref);
    }
}

void FluxSurfaceQuantity::DeallocateReferenceData(){
    if(referenceData==nullptr)
        return;
    
    for(len_t ir=0; ir<nr; ir++){
        delete [] referenceData[ir];
        gsl_spline_free(quantitySpline[ir]);
    }
    for(len_t ir=0; ir<=nr; ir++){
        delete [] referenceData_fr[ir];
        gsl_spline_free(quantitySpline_fr[ir]);
    }
    delete [] referenceData;
    delete [] referenceData_fr; 
}

const real_t *FluxSurfaceQuantity::GetData(len_t ir, fluxGridType fluxGridType) const {
    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        return quantityData_fr[ir];
    else
        return quantityData[ir];
}



const real_t FluxSurfaceQuantity::evaluateAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType) const {
    VerifyTheta(&theta);
    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        return gsl_spline_eval(quantitySpline_fr[ir], theta, gsl_acc);
    else
        return gsl_spline_eval(quantitySpline[ir], theta, gsl_acc);
}

// Shift into the interval for which there is magnetic field data,
// i.e. into theta_ref_min < theta < theta_ref_max
void FluxSurfaceQuantity::VerifyTheta(real_t *theta) const {
    if ( (*theta>=theta_ref_min) && (*theta<=theta_ref_max) )
        return;

    if (*theta < 0){
        // number of factors of 2pi to add to theta
        // to end up in [0,2pi]
        real_t n2Pi;
        std::modf(*theta/(-2*M_PI),&n2Pi);  
        *theta += 2*M_PI * (n2Pi+1);
    }

    if( *theta >= theta_ref_max )
        throw FVMException("Reference poloidal angle grid does not span a full orbit.");



}