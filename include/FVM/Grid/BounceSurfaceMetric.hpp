#ifndef _DREAM_FVM_BOUNCE_SURFACE_METRIC_HPP
#define _DREAM_FVM_BOUNCE_SURFACE_METRIC_HPP

namespace DREAM::FVM { class BounceSurfaceMetric; }

#include "FVM/Grid/BounceSurfaceQuantity.hpp"


namespace DREAM::FVM {
    class BounceSurfaceMetric : public BounceSurfaceQuantity {
    protected:
        FluxSurfaceQuantity *B;
        FluxSurfaceAverager *fluxSurfaceAverager;
        virtual void InterpolateToBounceGrid(real_t ***&bounceData, fluxGridType fluxGridType) override;
        
        void InterpolateToFluxGrid(
            real_t ***&bounceData, fluxGridType fluxGridType,
            len_t ntheta_interp_passing, const real_t *theta_passing);
        virtual void DeleteData(real_t ***&data, bool **isTrapped, len_t nr, len_t np1, len_t np2) override;
    public:
        BounceSurfaceMetric(Grid *grid, FluxSurfaceQuantity *Jacobian, FluxSurfaceQuantity *B, FluxSurfaceAverager *FSA);
        virtual ~BounceSurfaceMetric();

        void SetDataForPassing(len_t ntheta_interp_passing, const real_t *theta_passing);

        virtual const real_t *GetData(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const override;        

        const real_t evaluateAtTheta(len_t ir, len_t i, len_t j, real_t theta, fluxGridType) const;
        const real_t evaluateAtTheta(len_t ir, len_t i, len_t j, real_t theta, real_t cosTheta, real_t sinTheta, fluxGridType) const;
        virtual const real_t evaluateAtTheta(len_t, real_t, fluxGridType) const override 
            {FVMException("evaluateAtTheta(ir) invalid for Metric; use (ir,i,j)."); return 0;}

    };
}

#endif/*_DREAM_FVM_BOUNCE_SURFACE_METRIC_HPP*/

