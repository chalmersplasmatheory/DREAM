#ifndef _DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>
#include <gsl/gsl_spline.h>

namespace DREAM::FVM {
    class AnalyticBRadialGridGenerator : public RadialGridGenerator {
    private:
        // number of theta points that angle-dependent quantities are
        // evaluated on, and bounce/flux surface averaged over
        const len_t ntheta = 20;  
        const real_t *theta;
        const real_t *weightsTheta;
        real_t rMin, rMax, R0;
        len_t nrProfiles;
        real_t *rProfilesProvided, *GsProvided, *psisProvided, 
            *kappasProvided, *deltasProvided, *DeltasProvided;
        real_t *G, *psi, *kappa, *delta, *Delta,
            *GPrime, *psiPrime, *kappaPrime, *deltaPrime, *DeltaPrime;
        real_t *G_f, *psi_f, *kappa_f, *delta_f, *Delta_f,
            *GPrime_f, *psiPrime_f, *kappaPrime_f, *deltaPrime_f, *DeltaPrime_f;
        
        real_t *nabla_r2 ,*nabla_r2_f, *R, *R_f;
        
        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;
        real_t diffFunc(real_t r, std::function<real_t(real_t)> F); // = dF/dr at r

        static real_t effectivePassingFractionIntegrand(real_t, void*);
        real_t EvaluateBounceSurfaceIntegral(RadialGrid*, const MomentumGrid*, len_t, len_t, len_t, len_t, std::function<real_t(real_t,real_t)>);
        real_t EvaluateFluxSurfaceIntegral(RadialGrid*, len_t, bool, std::function<real_t(real_t)>);
        real_t EvaluateFluxSurfaceIntegral(RadialGrid*, len_t, bool, real_t*);
        void InterpolateInputProfileToGrid(real_t*, real_t*, real_t*,real_t*, real_t*, real_t*,real_t*);
        gsl_spline *spline_x;
        gsl_interp_accel *gsl_acc;
    public:
        AnalyticBRadialGridGenerator(const len_t nr,  real_t r0, 
             real_t ra,  real_t R0,
             real_t *, len_t , real_t *Gs, real_t *psi_p0s,
             real_t *kappas, real_t *deltas, real_t *Deltas);
        ~AnalyticBRadialGridGenerator();

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void RebuildJacobians(RadialGrid*, MomentumGrid**) override;
        virtual void RebuildFSAvgQuantities(RadialGrid*, MomentumGrid**) override;
        real_t FluxSurfaceAverageQuantity(RadialGrid *rGrid, len_t ir, bool rFluxGrid, std::function<real_t(real_t)> F) override
            {return EvaluateFluxSurfaceIntegral(rGrid, ir, rFluxGrid, F) / EvaluateFluxSurfaceIntegral(rGrid,ir,rFluxGrid, [](real_t){return 1;});} 
        real_t FluxSurfaceAverageQuantity(RadialGrid *rGrid, len_t ir, bool rFluxGrid, real_t *F) 
            {return EvaluateFluxSurfaceIntegral(rGrid, ir, rFluxGrid, F) / EvaluateFluxSurfaceIntegral(rGrid,ir,rFluxGrid, [](real_t){return 1;});} 
        virtual real_t BounceAverageQuantity(RadialGrid*, const MomentumGrid*, len_t, len_t, len_t, len_t, std::function<real_t(real_t,real_t)>) override;

        
    };
}

#endif/*_DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP*/
