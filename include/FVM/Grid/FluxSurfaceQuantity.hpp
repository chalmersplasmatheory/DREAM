#ifndef _DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP
#define _DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP

namespace DREAM::FVM { class FluxSurfaceQuantity; }

#include "FVM/Grid/RadialGrid.hpp"
#include "gsl/gsl_spline.h"


namespace DREAM::FVM {
    class FluxSurfaceQuantity {
    
    private:
        real_t 
            **quantityData     = nullptr,
            **quantityData_fr  = nullptr,
            **referenceData    = nullptr,
            **referenceData_fr = nullptr;

        RadialGrid *rGrid; 
        len_t nr;        

        gsl_spline  
            **quantitySpline = nullptr,
            **quantitySpline_fr = nullptr;
        

        real_t *theta_ref = nullptr;
        real_t theta_ref_min;
        real_t theta_ref_max;
        len_t ntheta_ref;

        //real_t *theta;
        //len_t ntheta;


        gsl_interp_accel *gsl_acc;
        const gsl_interp_type *interpolationMethod;



        void DeallocateReferenceData();
        void DeallocateInterpolatedData();
        void VerifyTheta(real_t *theta) const;
    public:
        FluxSurfaceQuantity(RadialGrid *rGrid, const gsl_interp_type *interpType);

        ~FluxSurfaceQuantity();

        
        void SetData(real_t *theta);
        const real_t *GetData(len_t ir, fluxGridType) const;
        const real_t evaluateAtTheta(len_t ir, real_t theta, fluxGridType) const;
        
        void Initialize(real_t **referenceData, real_t **referenceData_fr, real_t *theta_ref, len_t ntheta_ref);
        void InterpolateMagneticDataToTheta(real_t *theta, len_t ntheta_interp);

        gsl_spline *const* GetInterpolator() const {return quantitySpline;}
        gsl_spline *const* GetInterpolator_fr() const {return quantitySpline_fr;}
        real_t *const* GetData() const {return quantityData;}
        real_t *const* GetData_fr() const {return quantityData_fr;}
    };


}


#endif/*_DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP*/

