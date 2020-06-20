#ifndef _DREAM_FVM_BOUNCE_AVERAGER_HPP
#define _DREAM_FVM_BOUNCE_AVERAGER_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include <functional>
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"

namespace DREAM::FVM {
    class BounceAverager {
    private:
        // Pointer to the grid which owns this BounceAverager.
        Grid *grid;

        len_t nr;
        len_t *np1, *np2;

        bool geometryIsSymmetric;
        bool gridModePXI;
        /**
         * Is true if BounceAverager is constructed with 
         * QUAD_ADAPTIVE. Overrides a bunch of stuff and evaluates
         * flux surface averages with an adaptive quadrature.
         */ 
        bool integrateTrappedAdaptive = false;
        bool integratePassingAdaptive;
        // Pointer to the FluxSurfaceAverager owned by grid->radialgrid.
        FluxSurfaceAverager *fluxSurfaceAverager;

        // poloidal grid and quadrature weights for passing-particle grid
        // from fluxSurfaceAverager.
        len_t ntheta_interp;
        const real_t 
            *theta,
            *weights;

        // Reference quadrature grid points and weights on x in [0,1].
        len_t ntheta_interp_trapped;
        real_t *quad_x_ref;
        real_t *quad_w_ref;
        //std::function<real_t(real_t,real_t,real_t)> quadratureWeightFunction;
        gsl_integration_fixed_workspace *gsl_w = nullptr;
        gsl_integration_workspace *gsl_adaptive = nullptr;

        // True if particle is on a trapped orbit.
        // Size nr(+1) x np1(+1) x np2(+1)
        bool    **isTrapped    = nullptr,    // on distribution grid 
                **isTrapped_fr = nullptr,    // on radial flux grid 
                **isTrapped_f1 = nullptr,    // on p1 flux grid 
                **isTrapped_f2 = nullptr;    // on p2 flux grid

        // Size NR+ x (NP1+ x NP2+).
        // If isTrapped, contains bounce point theta_b1 or theta_b2,
        // otherwise empty.
        real_t  **theta_b1    = nullptr, // on distribution grid 
                **theta_b1_fr = nullptr, // on radial flux grid 
                **theta_b1_f1 = nullptr, // on p1 flux grid
                **theta_b1_f2 = nullptr, // on p2 flux grid
                **theta_b2    = nullptr, // on distribution grid 
                **theta_b2_fr = nullptr, // on radial flux grid 
                **theta_b2_f1 = nullptr, // on p1 flux grid
                **theta_b2_f2 = nullptr; // on p2 flux grid

        // Size NR+ x (NP1+ x NP2+) x ntheta_interp. 
        // If isTrapped, contains theta grid between theta_b1 and theta_b2,
        // otherwise empty. (between 0 and theta_b2 if isUpDownSymmetric)
        real_t  ***theta_bounceGrid    = nullptr, // on distribution grid 
                ***theta_bounceGrid_fr = nullptr, // on radial flux grid 
                ***theta_bounceGrid_f1 = nullptr, // on p1 flux grid
                ***theta_bounceGrid_f2 = nullptr; // on p2 flux grid
        // If isTrapped, contains quadrature weights corresponding to theta_bounceGrid.
        real_t  ***weights_bounceGrid    = nullptr, // on distribution grid 
                ***weights_bounceGrid_fr = nullptr, // on radial flux grid 
                ***weights_bounceGrid_f1 = nullptr, // on p1 flux grid
                ***weights_bounceGrid_f2 = nullptr; // on p2 flux grid

        // If isTrapped, contains magnetic field evaluated on theta_bounceGrid.
        real_t  ***B_bounceGrid    = nullptr, // on distribution grid 
                ***B_bounceGrid_fr = nullptr, // on radial flux grid 
                ***B_bounceGrid_f1 = nullptr, // on p1 flux grid
                ***B_bounceGrid_f2 = nullptr; // on p2 flux grid

        // If isTrapped, contains R/R0 evaluated on theta_bounceGrid.
        real_t  ***ROverR0_bounceGrid    = nullptr, // on distribution grid 
                ***ROverR0_bounceGrid_fr = nullptr, // on radial flux grid 
                ***ROverR0_bounceGrid_f1 = nullptr, // on p1 flux grid
                ***ROverR0_bounceGrid_f2 = nullptr; // on p2 flux grid

        // If isTrapped, contains Jacobian (1/R0 times the spatial 3D jacobian 
        // for r-theta-phi coordinates) evaluated on theta_bounceGrid.
        real_t  ***NablaR2_bounceGrid    = nullptr, // on distribution grid 
                ***NablaR2_bounceGrid_fr = nullptr, // on radial flux grid 
                ***NablaR2_bounceGrid_f1 = nullptr, // on p1 flux grid
                ***NablaR2_bounceGrid_f2 = nullptr; // on p2 flux grid


        // Size NR+ x (NP1+ x NP2+) x ntheta_interp.
        // Contains the full metric "J*sqrt(g)" on entire grid + poloidal.
        real_t  ***metricSqrtG    = nullptr,
                ***metricSqrtG_fr = nullptr,
                ***metricSqrtG_f1 = nullptr,
                ***metricSqrtG_f2 = nullptr;

        // Size NR+ x (NP1+ x NP2+) .
        // Contains the bounce-integrated metric VPrime on entire grid.
        real_t  **Vp    = nullptr,
                **Vp_fr = nullptr,
                **Vp_f1 = nullptr,
                **Vp_f2 = nullptr;


        virtual real_t EvaluateBounceIntegral(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F);

        const real_t* GetB(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t* GetROverR0(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t* GetNablaR2(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t* GetTheta(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t* GetWeights(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t* GetMetric(len_t ir, len_t i, len_t j, fluxGridType fluxGridType) const;
        const real_t GetVp(len_t, len_t, len_t, fluxGridType) const;
        const real_t GetIsTrapped(len_t, len_t, len_t, fluxGridType) const;

        virtual void FindBouncePoints(len_t ir, real_t xi0, bool rFluxGrid, real_t *thetab_1, real_t *thetab_2);
        static real_t xiParticleFunction(real_t, void*);
        virtual void FindThetaBounceRoots(real_t *x_lo, real_t *x_up, real_t *root, gsl_function);

        void InitializeQuadrature(FluxSurfaceAverager::quadrature_method);
        bool SetIsTrapped(bool**&, real_t**&, real_t**&, fluxGridType);

        bool InitializeBounceIntegralQuantities();
        void InterpolateMagneticQuantitiesToBounceGrid(
            real_t ***&, real_t ***&, real_t**, real_t**, real_t ***&, 
            real_t ***&, real_t ***&, real_t ***&, fluxGridType);

        void SetMetricOnPassingGrid(real_t***&, bool**, fluxGridType);
        void SetVp(real_t**&, fluxGridType);

        void AllocateBounceIntegralQuantities();
        void AllocateBounceGridMagneticQuantities(
            real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid,  real_t ***&B_bounceGrid, 
            real_t ***&ROverR0_bounceGrid, real_t ***&NablaR2_bounceGrid, fluxGridType fluxGridType
        );
        void DeallocateBounceIntegralQuantities();
        void DeallocateBounceGridMagneticQuantities(            
            real_t ***&theta_bounceGrid, real_t ***&weights_bounceGrid,  real_t ***&B_bounceGrid, 
            real_t ***&ROverR0_bounceGrid, real_t ***&NablaR2_bounceGrid, fluxGridType fluxGridType
        );

        gsl_interp_accel *gsl_acc;

    public:
        BounceAverager(
            Grid*, FluxSurfaceAverager*, len_t ntheta_interp_trapped,
            enum OptionConstants::momentumgrid_type mgtype,
            FluxSurfaceAverager::quadrature_method q_method_trapped = FluxSurfaceAverager::QUAD_FIXED_LEGENDRE
        );
        
        ~BounceAverager();

        virtual real_t CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F);

        void Rebuild();

    };
}

#endif/*_DREAM_FVM_BOUNCE_AVERAGER_HPP*/
