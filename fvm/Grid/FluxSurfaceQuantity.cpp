/**
 * Implementation of the FluxSurfaceQuantity class which contains data and calculations
 * of poloidal-angle-dependent quantities used in flux surface averages. Contains
 * interpolators that can interpolate quantity to any theta, and methods to store 
 * on grid for fixed quadratures.
 */

#include "FVM/Grid/FluxSurfaceQuantity.hpp"

using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor
 */
FluxSurfaceQuantity::FluxSurfaceQuantity(
    RadialGrid *rGrid, std::function<real_t(len_t,real_t)> evalAtTheta, 
    std::function<real_t(len_t,real_t)> evalAtTheta_f, 
    const gsl_interp_type *interpType
) : rGrid(rGrid), evaluateFuncAtTheta(evalAtTheta), evaluateFuncAtTheta_f(evalAtTheta_f),  interpolationMethod(interpType) 
{ 
    gsl_acc = gsl_interp_accel_alloc();
}

/**
 * Destructor
 */
FluxSurfaceQuantity::~FluxSurfaceQuantity(){
    gsl_interp_accel_free(gsl_acc);
    DeallocateInterpolatedData();
}

/**
 * Interpolate (and store) data to theta grid.
 */
void FluxSurfaceQuantity::InterpolateMagneticDataToTheta(real_t *theta, len_t ntheta_interp){
    DeallocateInterpolatedData();
    nr = rGrid->GetNr();
    quantityData = new real_t*[nr];
    for(len_t ir = 0; ir<nr; ir++){
        quantityData[ir] = new real_t[ntheta_interp];
        for(len_t it=0; it<ntheta_interp; it++)
            quantityData[ir][it] = evaluateFuncAtTheta(ir,theta[it]);
    }
    quantityData_fr = new real_t*[nr+1];
    for(len_t ir = 0; ir<=nr; ir++){
        quantityData_fr[ir] = new real_t[ntheta_interp];
        for(len_t it=0; it<ntheta_interp; it++)
            quantityData_fr[ir][it] = evaluateFuncAtTheta_f(ir,theta[it]);
    }
}

/**
 * Deallocator
 */
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


/**
 * Get stored data (when using fixed quadratures) 
 */
const real_t *FluxSurfaceQuantity::GetData(len_t ir, fluxGridType fluxGridType) const {
    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        return quantityData_fr[ir];
    else
        return quantityData[ir];
}

/**
 * Evaluate quantity at any poloidal angle theta
 */
const real_t FluxSurfaceQuantity::evaluateAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType) const {
    if (fluxGridType == FLUXGRIDTYPE_RADIAL)
        return evaluateFuncAtTheta_f(ir,theta);
    else
        return evaluateFuncAtTheta(ir,theta);
}
