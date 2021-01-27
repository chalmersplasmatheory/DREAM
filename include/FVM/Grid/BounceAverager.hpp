#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/BounceSurfaceMetric.hpp"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_roots.h"

#ifndef _DREAM_FVM_BOUNCE_AVERAGER_HPP
#define _DREAM_FVM_BOUNCE_AVERAGER_HPP

namespace DREAM::FVM {
    class BounceAverager {
    private:
        // Pointer to the grid which owns this BounceAverager.
        Grid *grid;

        // Resolution parameters.
        len_t nr,
            *np1 = nullptr, 
            *np2 = nullptr;

        bool geometryIsSymmetric;

        // Pointer to the FluxSurfaceAverager owned by grid->radialgrid.
        FluxSurfaceAverager *fluxSurfaceAverager;

        // poloidal grid and quadrature weights for passing-particle grid
        // from fluxSurfaceAverager.
        len_t ntheta_interp_passing;
        const real_t 
            *theta_passing,
            *weights_passing;

        // Reference quadrature grid points and weights on x in [0,1].
        len_t ntheta_interp_trapped;
        real_t *theta_trapped_ref;
        real_t *weights_trapped_ref;

        // true if evaluate bounce integral with adaptive quadrature
        bool integrateTrappedAdaptive = false; //...for trapped orbits
        bool integratePassingAdaptive = false; //...for passing orbits

        gsl_integration_fixed_workspace *gsl_w = nullptr;
        gsl_integration_workspace *gsl_adaptive = nullptr;
        gsl_root_fsolver *gsl_fsolver = nullptr;
        gsl_integration_qaws_table *qaws_table;
        int QAG_KEY = GSL_INTEG_GAUSS41;
        
        BounceSurfaceQuantity *BOverBmin;
        BounceSurfaceQuantity *ROverR0;
        BounceSurfaceQuantity *NablaR2;
        BounceSurfaceMetric   *Metric;
        
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

        real_t EvaluateBounceIntegralOverP2(len_t ir, len_t i, len_t j, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *F_list=nullptr);
        void InitializeQuadrature(FluxSurfaceAverager::quadrature_method);
        bool SetIsTrapped(bool**&, real_t**&, real_t**&, fluxGridType);

        bool InitializeBounceIntegralQuantities();
        void SetVp(real_t**&, fluxGridType);
        void SetVpsPXi(real_t**&,real_t**&,real_t**&);

        void AllocateBounceIntegralQuantities();

        real_t GetXi0(len_t ir, len_t i, len_t j, fluxGridType);
        real_t GetVp(len_t ir, len_t i, len_t j, fluxGridType);
        real_t GetBmin(len_t ir, fluxGridType);
        real_t GetBmax(len_t ir, fluxGridType);
        
        void UpdateGridResolution();

        const real_t realeps = std::numeric_limits<real_t>::epsilon();
    public:
        BounceAverager(
            Grid*, FluxSurfaceAverager*, len_t ntheta_interp_trapped,
            FluxSurfaceAverager::quadrature_method q_method_trapped = FluxSurfaceAverager::QUAD_FIXED_CHEBYSHEV
        );
        ~BounceAverager();

        real_t CalculateBounceAverage(len_t ir, len_t i, len_t j, fluxGridType fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par, const int_t *F_list=nullptr);
        void Rebuild();

        BounceSurfaceQuantity *GetBOverBmin(){return BOverBmin;}
        BounceSurfaceQuantity *GetROverR0(){return ROverR0;}
        BounceSurfaceQuantity *GetNablaR2(){return NablaR2;}
        BounceSurfaceMetric   *GetMetric(){return Metric;}
        FluxSurfaceAverager *GetFluxSurfaceAverager(){return fluxSurfaceAverager;}        
    };
}

#endif/*_DREAM_FVM_BOUNCE_AVERAGER_HPP*/
