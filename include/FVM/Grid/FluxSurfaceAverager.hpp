#ifndef _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP
#define _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP

namespace DREAM::FVM { class FluxSurfaceAverager; }

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/FluxSurfaceQuantity.hpp"
#include <functional>
#include "gsl/gsl_integration.h"




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
        bool integrateAdaptive;

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

        gsl_integration_fixed_workspace *gsl_w;
        gsl_integration_workspace *gsl_adaptive;
//        real_t *quadratureWeightFunction;
        len_t ntheta_interp; // number of poloidal grid points
        real_t  
            *theta   = nullptr, // poloidal grid points
            *weights = nullptr, // corresponding quadrature weights
            theta_max;


        void InitializeQuadrature(quadrature_method);
        void DeallocateQuadrature();
        
        static real_t FluxSurfaceIntegralFunction(real_t x, void *p);


    public:
        FluxSurfaceAverager(
            RadialGrid*, bool geometryIsSymmetric = false, len_t ntheta_interp = 10,
            interp_method = INTERP_LINEAR, quadrature_method = QUAD_FIXED_LEGENDRE
        );
        ~FluxSurfaceAverager();

        void Rebuild();

        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType, std::function<real_t(real_t,real_t,real_t)>);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType, std::function<real_t(real_t,real_t,real_t)>);

        const len_t GetNTheta() const
            {return ntheta_interp;}    
        const real_t *GetTheta() const
            {return theta;}        
        const real_t *GetWeights() const
            {return weights;}        

        real_t evaluateXiAtTheta(len_t ir, real_t xi0, real_t theta, fluxGridType fluxGridType);

        static real_t evaluateXiAtB(real_t xi0, real_t BOverBmin){
            real_t sgnXi0 = (xi0>0) - (xi0<0);
            return sgnXi0 * sqrt(1- (1-xi0*xi0)*BOverBmin );
        }

        void SetReferenceMagneticFieldData(
            len_t ntheta_ref, real_t *theta_ref,
            real_t **B_ref, real_t **B_ref_f,
            real_t **Jacobian_ref, real_t **Jacobian_ref_f,
            real_t **ROverR0_ref ,real_t **ROverR0_ref_f, 
            real_t **NablaR2_ref, real_t **NablaR2_ref_f
        );

        bool isGeometrySymmetric(){return geometryIsSymmetric;}
        bool isIntegrationAdaptive(){return integrateAdaptive;}

    };
}

#endif/*_DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP*/

