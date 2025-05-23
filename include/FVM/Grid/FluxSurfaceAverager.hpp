#include "FVM/Grid/RadialGrid.hpp"

#ifndef _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP
#define _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP

namespace DREAM::FVM { class FluxSurfaceAverager; }

#include "FVM/Grid/FluxSurfaceQuantity.hpp"
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

        // parameters for BounceIntegralFunction
        struct BounceIntegralParams {
            len_t ir; real_t xi0; real_t theta_b1; real_t theta_b2; fluxGridType fgType; 
            real_t Bmin; real_t(*F_ref)(real_t,real_t,real_t,real_t,void*); real_t(*F_eval)(real_t,real_t,real_t,real_t,void*); void *F_ref_par; int_t *Flist_eval; 
            FluxSurfaceAverager *fsAvg; bool integrateQAWS;
        };
        static real_t BA_FUNC_PASSING(real_t xiOverXi0, real_t BOverBmin, real_t ROverR0, real_t NablaR2, void* par){
            BounceIntegralParams *params = (BounceIntegralParams*)par;
            return params->F_ref(xiOverXi0, BOverBmin, ROverR0, NablaR2, params->F_ref_par);
        }
        static real_t BA_FUNC_TRAPPED(real_t xiOverXi0, real_t BOverBmin, real_t ROverR0, real_t NablaR2, void* par){
            BounceIntegralParams *params = (BounceIntegralParams*)par;
            return params->F_ref( xiOverXi0, BOverBmin, ROverR0, NablaR2, params->F_ref_par) 
                 + params->F_ref(-xiOverXi0, BOverBmin, ROverR0, NablaR2, params->F_ref_par);
        }

    private:
        // Pointer to the RadialGrid which owns 
        // this FluxSurfaceAverager, and its generator.
        RadialGrid *rGrid;
        RadialGridGenerator *gridGenerator;

        bool geometryIsSymmetric;

        /**
         * Is true if FluSurfaceAverager is constructed with 
         * QUAD_ADAPTIVE. Overrides a bunch of stuff and evaluates
         * flux surface averages with an adaptive quadrature.
         */ 
        bool integrateAdaptive = false;

        // Number of radial grid points on distribution grid
        len_t nr;

        FluxSurfaceQuantity 
            *BOverBmin = nullptr,
            *Jacobian = nullptr,
            *ROverR0 = nullptr,
            *NablaR2 = nullptr;
        
        gsl_integration_fixed_workspace *gsl_w = nullptr;
        gsl_integration_workspace *gsl_adaptive;
        gsl_integration_workspace *gsl_adaptive_outer;
        gsl_root_fsolver *gsl_fsolver;
        gsl_integration_qaws_table *qaws_table;
        int QAG_KEY = GSL_INTEG_GAUSS41;

        len_t ntheta_interp; // number of poloidal grid points
        real_t  
            *theta   = nullptr, // poloidal grid points
            *weights = nullptr, // corresponding quadrature weights
            theta_max;

        // poloidal angles of minimum and maximum magnetic field strength.
        real_t 
            *theta_Bmin = nullptr,
            *theta_Bmin_f,
            *theta_Bmax,
            *theta_Bmax_f;

        void InitializeQuadrature(quadrature_method);
        void DeallocateQuadrature();

        void InitializeReferenceData(
            real_t *theta_Bmin, real_t *theta_Bmin_f,
            real_t *theta_Bmax, real_t *theta_Bmax_f
        );
        void DeallocateReferenceData();
        
        static real_t FluxSurfaceIntegralFunction(real_t x, void *p);

        real_t GetVpVol(len_t ir, fluxGridType);

    static real_t evaluatePXiBounceIntegralAtXi(real_t,void*);

    public:
        FluxSurfaceAverager(
            RadialGrid*, RadialGridGenerator*, bool geometryIsSymmetric = false, len_t ntheta_interp = 10,
            interp_method im = INTERP_LINEAR, quadrature_method qm = QUAD_FIXED_LEGENDRE
        );
        ~FluxSurfaceAverager();

        void Rebuild();

        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t EvaluatePXiBounceIntegralAtP(len_t ir, real_t xi0, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t CalculatePXiBounceAverageAtP(len_t ir, real_t xi0, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);

        real_t EvaluateCellAveragedBounceIntegralOverP2(len_t ir, real_t xi_f1, real_t xi_f2, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        bool shouldCellAverageBounceIntegral(len_t ir, real_t xi_lower, real_t xi_upper, fluxGridType);

        static real_t BounceIntegralFunction(real_t theta, void *par);

        const len_t GetNTheta() const
            {return ntheta_interp;}    
        const real_t *GetTheta() const
            {return theta;}        
        const real_t *GetWeights() const
            {return weights;}        
        const real_t GetThetaMax() const    
            {return theta_max;}

        void SetReferenceMagneticFieldData(
            real_t *theta_Bmin, real_t *theta_Bmin_f, // poloidal angle of B=Bmin
            real_t *theta_Bmax, real_t *theta_Bmax_f  // poloidal angle of B=Bmax
        );

        FluxSurfaceQuantity *GetBOverBmin(){return BOverBmin;}
        FluxSurfaceQuantity *GetJacobian(){return Jacobian;}
        FluxSurfaceQuantity *GetROverR0(){return ROverR0;}
        FluxSurfaceQuantity *GetNablaR2(){return NablaR2;}
        real_t BAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->BAtTheta_f(ir,theta);
            else
                return gridGenerator->BAtTheta(ir,theta);            
        }
        real_t JacobianAtTheta(len_t ir, real_t theta, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->JacobianAtTheta_f(ir,theta);
            else
                return gridGenerator->JacobianAtTheta(ir,theta);            
        }
        real_t ROverR0AtTheta(len_t ir, real_t theta, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->ROverR0AtTheta_f(ir,theta);
            else
                return gridGenerator->ROverR0AtTheta(ir,theta);            
        }
        real_t NablaR2AtTheta(len_t ir, real_t theta, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->NablaR2AtTheta_f(ir,theta);
            else
                return gridGenerator->NablaR2AtTheta(ir,theta);            
        }
        void GeometricQuantitiesAtTheta(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                gridGenerator->EvaluateGeometricQuantities_fr(ir, theta, B, Jacobian, ROverR0, NablaR2);
            else 
                gridGenerator->EvaluateGeometricQuantities(ir, theta, B, Jacobian, ROverR0, NablaR2);
        }

        real_t GetBmin(len_t ir, fluxGridType, real_t *theta_Bmin = nullptr);
        real_t GetBmax(len_t ir, fluxGridType, real_t *theta_Bmax = nullptr);


        bool isGeometrySymmetric(){return geometryIsSymmetric;}
        bool isIntegrationAdaptive(){return integrateAdaptive;}


        static void FindThetas(real_t theta_Bmin, real_t theta_Bmax, real_t *theta1, real_t *theta2, gsl_function, gsl_root_fsolver*, bool isSymmetric=false);
        static void FindRoot(real_t *x_lo, real_t *x_up, real_t *root, gsl_function, gsl_root_fsolver*, real_t epsrel=1e-7, real_t epsabs=1e-7, len_t max_iter=50);
        static void FindBouncePoints(len_t ir, real_t Bmin, real_t theta_Bmin, real_t theta_Bmax, FluxSurfaceAverager*, real_t xi0, fluxGridType, real_t *thetab_1, real_t *thetab_2, gsl_root_fsolver*, bool isSymmetric=false);
        static real_t xiParticleFunction(real_t, void*);

        static real_t AssembleBAFunc(
            real_t xiOverXi0,real_t BOverBmin, real_t ROverR0, 
            real_t NablaR2, const int_t *Flist
        );
        static real_t AssembleFSAFunc(
            real_t BOverBmin, real_t ROverR0, real_t NablaR2, const int_t *Flist
        );
        real_t EvaluateAvalancheDeltaHat(len_t ir, real_t p, real_t xi_l, real_t xi_u, real_t Vp, real_t VpVol, int_t RESign = 1);

		void PrintBOfTheta(const len_t ir, const len_t N=1000, enum fluxGridType fgt=FLUXGRIDTYPE_DISTRIBUTION);
    };
}

#endif/*_DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP*/

