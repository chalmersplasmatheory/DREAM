#include "FVM/Grid/RadialGrid.hpp"

#ifndef _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP
#define _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP

namespace DREAM::FVM { class FluxSurfaceAverager; }

#include "FVM/Grid/FluxSurfaceQuantity.hpp"
#include <functional>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_roots.h"




namespace DREAM::FVM {
    class FluxSurfaceAverager {

    public:
        /**
         * The interpolation method used when
         * interpolating from reference magnetic
         * data from RadialGridGenerator to grid
         * points requested by quadrature rule.
         */ 
        enum interp_method {
            INTERP_LINEAR, // Linear interpolation
            INTERP_STEFFEN // Extremum-preserving spline
        };

        /**
         * Which quadrature rule is used for 
         * evaluating flux surface integrals.
         */
        enum quadrature_method {
            QUAD_FIXED_LEGENDRE,
            QUAD_FIXED_CHEBYSHEV,
            QUAD_ADAPTIVE
        };

    private:
        // Pointer to the RadialGrid which owns 
        // this FluxSurfaceAverager.
        RadialGrid *rGrid;


        /**
         * True if geometry is up-down symmetric, and 
         * theta_ref[0] = 0, theta_ref[ntheta_ref-1] = pi
         * and B[ir][0] = Bmin[ir] (i.e. 0 is outboard side).
         */
        bool geometryIsSymmetric;

        /**
         * Is true if FluSurfaceAverager is constructed with 
         * QUAD_ADAPTIVE. Overrides a bunch of stuff and evaluates
         * flux surface averages with an adaptive quadrature.
         */ 
        bool integrateAdaptive = false;

        // Number of radial grid points on distribution grid
        len_t nr;


        // Number of reference poloidal angle points
        len_t ntheta_ref;

        /**
         * Poloidal angle grid points on reference
         * grid used by RadialGridGenerator
         * (same on all radii. Size: ntheta_ref)
         */
        real_t *theta_ref = nullptr;

        FluxSurfaceQuantity 
            *B = nullptr,
            *Jacobian = nullptr,
            *ROverR0 = nullptr,
            *NablaR2 = nullptr;

        gsl_integration_fixed_workspace *gsl_w = nullptr;
        gsl_integration_workspace *gsl_adaptive;
        gsl_root_fsolver *gsl_fsolver;
//        real_t *quadratureWeightFunction;
        len_t ntheta_interp; // number of poloidal grid points
        real_t  
            *theta   = nullptr, // poloidal grid points
            *weights = nullptr, // corresponding quadrature weights
            theta_max;


        void InitializeQuadrature(quadrature_method);
        void DeallocateQuadrature();
        
        static real_t FluxSurfaceIntegralFunction(real_t x, void *p);

        real_t GetBmin(len_t ir, fluxGridType);
        real_t GetBmax(len_t ir, fluxGridType);
        real_t GetVpVol(len_t ir, fluxGridType);

        static void FindThetaBounceRoots(real_t *x_lo, real_t *x_up, real_t *root, gsl_function, gsl_root_fsolver*);

    public:
        FluxSurfaceAverager(
            RadialGrid*, bool geometryIsSymmetric = false, len_t ntheta_interp = 10,
            interp_method im = INTERP_LINEAR, quadrature_method qm = QUAD_FIXED_LEGENDRE
        );
        ~FluxSurfaceAverager();

        void Rebuild();

        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType, std::function<real_t(real_t,real_t,real_t)>);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType, std::function<real_t(real_t,real_t,real_t)>);
        real_t EvaluatePXiBounceIntegralAtP(len_t ir, real_t p, real_t xi0, fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F);
        real_t CalculatePXiBounceAverageAtP(len_t ir, real_t p, real_t xi0, fluxGridType fluxGridType, std::function<real_t(real_t,real_t,real_t,real_t)> F);

        const len_t GetNTheta() const
            {return ntheta_interp;}    
        const real_t *GetTheta() const
            {return theta;}        
        const real_t *GetWeights() const
            {return weights;}        
        const real_t GetThetaMax() const    
            {return theta_max;}

        void SetReferenceMagneticFieldData(
            len_t ntheta_ref, real_t *theta_ref,
            real_t **B_ref, real_t **B_ref_f,
            real_t **Jacobian_ref, real_t **Jacobian_ref_f,
            real_t **ROverR0_ref ,real_t **ROverR0_ref_f, 
            real_t **NablaR2_ref, real_t **NablaR2_ref_f
        );

        FluxSurfaceQuantity *GetB(){return B;}
        FluxSurfaceQuantity *GetJacobian(){return Jacobian;}
        FluxSurfaceQuantity *GetROverR0(){return ROverR0;}
        FluxSurfaceQuantity *GetNablaR2(){return NablaR2;}
        

        bool isGeometrySymmetric(){return geometryIsSymmetric;}
        bool isIntegrationAdaptive(){return integrateAdaptive;}



        static void FindBouncePoints(len_t ir, real_t Bmin, const FluxSurfaceQuantity *B, real_t xi0, fluxGridType, real_t *thetab_1, real_t *thetab_2, gsl_root_fsolver*, bool geometryIsSymmetric = false);
        static real_t xiParticleFunction(real_t, void*);

    };
}

#endif/*_DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP*/

