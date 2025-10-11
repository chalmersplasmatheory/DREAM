#include "FVM/Grid/Stellarator/RadialGrid.hpp"

#ifndef _DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP // TODO
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
        /** TODO: Take back if BA
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
        }*/

    private:
        // Pointer to the RadialGrid which owns 
        // this FluxSurfaceAverager, and its generator.
        RadialGrid *rGrid;
        RadialGridGenerator *gridGenerator;

        len_t nfp = 0;

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
            *BdotGradphi = nullptr,
            *gttOverJ2 = nullptr,
            *gtpOverJ2 = nullptr;
        
        gsl_integration_fixed_workspace *gsl_w_theta = nullptr;
        gsl_integration_fixed_workspace *gsl_w_phi    = nullptr;
        gsl_integration_workspace *gsl_adaptive_theta;
        gsl_integration_workspace *gsl_adaptive_phi;
        // TODO: Take back if BA
        //gsl_integration_workspace *gsl_adaptive;
        //gsl_integration_workspace *gsl_adaptive_outer;
        //gsl_root_fsolver *gsl_fsolver;
        //gsl_integration_qaws_table *qaws_table;
        int QAG_KEY = GSL_INTEG_GAUSS41;

        len_t ntheta_interp, nphi_interp; // number of poloidal and toroidal grid points
        real_t  
            *theta   = nullptr, // poloidal grid points
            *phi     = nullptr, // toroidal grid points
            *weights_theta = nullptr, // corresponding quadrature weights for theta
            *weights_phi   = nullptr, // corresponding quadrature weights for phi
            theta_max = 2 * M_PI, phi_max = 2 * M_PI;

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
        
        static real_t FluxSurfaceIntegralFunctionThetaPhi(real_t,real_t, void *p);
        static real_t FluxSurfaceIntegralFunctionPhi(real_t, void *p);

        real_t GetVpVol(len_t ir, fluxGridType);

    public:
        FluxSurfaceAverager(
            RadialGrid*, RadialGridGenerator*, len_t nfp = 0, len_t ntheta_interp = 64, len_t nphi_interp = 64, // TODO: Change?
            interp_method im = INTERP_LINEAR, quadrature_method qm = QUAD_FIXED_LEGENDRE
        );
        ~FluxSurfaceAverager();

        void Rebuild();

        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        
        const len_t GetNTheta() const
            {return ntheta_interp;}    
        const real_t *GetTheta() const
            {return theta;}        
        const real_t *GetWeights() const
            {return weights;}        
        const real_t GetThetaMax() const    
            {return theta_max;}      
        const real_t GetPhiMax() const    
            {return phi_max;}

        void SetReferenceMagneticFieldData(
            real_t *theta_Bmin, real_t *theta_Bmin_f, // poloidal angle of B=Bmin
            real_t *theta_Bmax, real_t *theta_Bmax_f  // poloidal angle of B=Bmax
        );

        FluxSurfaceQuantity *GetBOverBmin(){return BOverBmin;}
        FluxSurfaceQuantity *GetJacobian(){return Jacobian;}
        FluxSurfaceQuantity *GetBdotGradphi(){return BdotGradphi;}
        FluxSurfaceQuantity *GetgttOverJ2(){return gttOverJ2;}
        FluxSurfaceQuantity *GetgtpOverJ2(){return gtpOverJ2;}
        real_t BAtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->BAtThetaPhi_f(ir,theta,phi);
            else
                return gridGenerator->BAtThetaPhi(ir,theta,phi);            
        }
        real_t JacobianAtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->JacobianAtThetaPhi_f(ir,theta,phi);
            else
                return gridGenerator->JacobianAtThetaPhi(ir,theta,phi);            
        }
        real_t BdotGradphiAtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->BdotGradphiAtThetaPhi_f(ir,theta);
            else
                return gridGenerator->BdotGradphiAtThetaPhi(ir,theta);            
        }
        real_t gttOverJ2AtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->gttAtThetaPhi_f(ir,theta);
            else
                return gridGenerator->gttAtThetaPhi(ir,theta);            
        }
        real_t gtpOverJ2AtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->gtpAtThetaPhi_f(ir,theta);
            else
                return gridGenerator->gtpAtThetaPhi(ir,theta);            
        }
        void GeometricQuantitiesAtThetaPhi(const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                gridGenerator->EvaluateGeometricQuantities(ir, theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2);
            else 
                gridGenerator->EvaluateGeometricQuantities_fr(ir, theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2);
        }

        real_t GetBmin(len_t ir, fluxGridType, real_t *theta_Bmin = nullptr);
        real_t GetBmax(len_t ir, fluxGridType, real_t *theta_Bmax = nullptr);


        len_t getNFP(){return nfp;}
        bool isIntegrationAdaptive(){return integrateAdaptive;}


        static real_t AssembleFSAFunc(
            real_t BOverBmin, real_t BdotGradphi, real_t gttOverJ2, real_t gtpOverJ2, const int_t *Flist
        );
        
		void PrintBOfThetaPhi(const len_t ir, const len_t N=100, enum fluxGridType fgt=FLUXGRIDTYPE_DISTRIBUTION);
    };
}

#endif/*_DREAM_FVM_FLUX_SURFACE_AVERAGER_HPP*/

