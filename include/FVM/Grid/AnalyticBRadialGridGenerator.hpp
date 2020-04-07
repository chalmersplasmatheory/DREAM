#ifndef _DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include <functional>

namespace DREAM::FVM {
    class AnalyticBRadialGridGenerator : public RadialGridGenerator {
    private:
        // number of theta points that angle-dependent quantities are
        // evaluated on, and bounce/flux surface averaged over
        const len_t ntheta = 100;  
        real_t theta;
        real_t rMin=0, rMax=1, R0=3;
        std::function<real_t(real_t)> 
            Btor_G       = [](real_t){return 3.0;}, // default G=3 Tm so that a default toroidal field strength B_phi = G/R0 = 1 T on axis.
            Psi_p0       = [](real_t){return 0.0;}, // default pure 1/R toroidal field
            kappa_Elong  = [](real_t){return 1.0;}, // default circular plasma
            delta_triang = [](real_t){return 0.0;}, // default no triangularity
            Delta_shafr  = [](real_t){return 0.0;}; // default no Shafranov shift
        
        //real_t* JacobianJ;
        //real_t JacobianJ_f;
        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;
        real_t diffFunc(real_t r, std::function<real_t(real_t)> F); // = dF/dr at r

        real_t EvaluateBounceSurfaceIntegral(RadialGrid*, const MomentumGrid*, len_t, len_t, len_t, len_t, std::function<real_t(real_t,real_t)>);
        real_t EvaluateFluxSurfaceIntegral(RadialGrid*, len_t  , len_t , std::function<real_t(real_t)>);

    public:
        AnalyticBRadialGridGenerator( len_t nr,  real_t r0, 
             real_t ra,  real_t R0,
            std::function<real_t(real_t)> G,  std::function<real_t(real_t)> Psi_p0, 
            std::function<real_t(real_t)> kappa, std::function<real_t(real_t)> delta, 
            std::function<real_t(real_t)> Delta);
        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void RebuildJacobians(RadialGrid*, MomentumGrid**) override;
        virtual void RebuildFSAvgQuantities(RadialGrid*, MomentumGrid**) override;
        virtual real_t FluxSurfaceAverageQuantity(RadialGrid*, len_t, bool, std::function<real_t(real_t)>) override;
        virtual real_t BounceAverageQuantity(RadialGrid*, const MomentumGrid*, len_t, len_t, len_t, len_t, std::function<real_t(real_t,real_t)>) override;

    };
}

#endif/*_DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP*/
