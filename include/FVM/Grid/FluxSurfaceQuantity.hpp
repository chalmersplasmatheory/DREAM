#ifndef _DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP
#define _DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP

namespace DREAM::FVM { class FluxSurfaceQuantity; }

#include "FVM/Grid/RadialGrid.hpp"
#include "gsl/gsl_spline.h"
#include <algorithm>

namespace DREAM::FVM {
    class FluxSurfaceQuantity {
    
    private:
        real_t 
            **quantityData     = nullptr,
            **quantityData_fr  = nullptr;

        RadialGrid *rGrid; 
        len_t nr;        

        gsl_spline  
            **quantitySpline = nullptr,
            **quantitySpline_fr = nullptr;
        
        std::function<real_t(len_t,real_t)> evaluateFuncAtTheta;
        std::function<real_t(len_t,real_t)> evaluateFuncAtTheta_f;

        gsl_interp_accel *gsl_acc;
        const gsl_interp_type *interpolationMethod;

        void DeallocateInterpolatedData();
    public:
        FluxSurfaceQuantity(
            RadialGrid *rGrid, std::function<real_t(len_t,real_t)>, std::function<real_t(len_t,real_t)>, 
            const gsl_interp_type *interpType
            );

        ~FluxSurfaceQuantity();
        
        void SetData(real_t *theta);
        const real_t *GetData(len_t ir, fluxGridType) const;
        const real_t evaluateAtTheta(len_t ir, real_t theta, fluxGridType) const;
        
        void InterpolateMagneticDataToTheta(real_t *theta, len_t ntheta_interp);

        real_t *const* GetData() const {return quantityData;}
        real_t *const* GetData_fr() const {return quantityData_fr;}
    };


}


#endif/*_DREAM_FVM_FLUX_SURFACE_QUANTITY_HPP*/

