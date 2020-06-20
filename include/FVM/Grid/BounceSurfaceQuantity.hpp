#ifndef _DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP
#define _DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP

namespace DREAM::FVM { class BounceSurfaceQuantity; }

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/FluxSurfaceQuantity.hpp"
#include "gsl/gsl_spline.h"


namespace DREAM::FVM {
    class BounceSurfaceQuantity {
    
    private:

        real_t
            ***bounceData      = nullptr,
            ***bounceData_fr   = nullptr,
            ***bounceData_f1   = nullptr,
            ***bounceData_f2   = nullptr;

        bool 
            **isTrapped = nullptr, 
            **isTrapped_fr = nullptr, 
            **isTrapped_f1 = nullptr,
            **isTrapped_f2 = nullptr;


        FluxSurfaceQuantity *fluxSurfaceQuantity;

        Grid *grid; 
        len_t nr;        
        len_t *np1;
        len_t *np2;


        gsl_interp_accel *gsl_acc;
        void DeallocateBounceData();

        const bool IsTrapped(len_t ir, len_t i, len_t j, fluxGridType) const;
        const real_t *GetBounceData(len_t ir, len_t i, len_t j, fluxGridType) const;
    public:
        BounceSurfaceQuantity(Grid *grid, FluxSurfaceQuantity *fluxSurfaceQuantity);
        ~BounceSurfaceQuantity();

        void SetData(real_t *theta);
        const real_t *GetData(len_t ir, len_t i, len_t j, fluxGridType) const;
        const real_t evaluateAtTheta(len_t ir, real_t theta, fluxGridType) const;
        
        void Initialize(bool **isTrapped, bool **isTrapped_fr, 
            bool **isTrapped_f1, bool **isTrapped_f2 );
            
        void InterpolateMagneticDataToBounceGrids(
            len_t ntheta_interp_trapped, real_t ***bounceTheta, real_t ***bounceTheta_fr,
            real_t ***bounceTheta_f1, real_t ***bounceTheta_f2);

    };


}


#endif/*_DREAM_FVM_BOUNCE_SURFACE_QUANTITY_HPP*/

