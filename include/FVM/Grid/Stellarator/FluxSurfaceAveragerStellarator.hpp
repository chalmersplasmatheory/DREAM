#include "FVM/Grid/Stellarator/RadialGridStellarator.hpp"

#ifndef _DREAM_FVM_FLUX_SURFACE_AVERAGER_STELLARATOR_HPP
#define _DREAM_FVM_FLUX_SURFACE_AVERAGER_STELLARATOR_HPP

namespace DREAM::FVM { class FluxSurfaceAveragerStellarator; }

#include "FVM/Grid/FluxSurfaceQuantity.hpp"
#include "FVM/Grid/FluxSurfaceAverager.hpp"
//#include "FVM/Grid/Stellarator/NumericStellaratorRadialGridGenerator.hpp"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_roots.h"


namespace DREAM::FVM {
    class FluxSurfaceAveragerStellarator : public FluxSurfaceAverager {

    public:

        // parameters for BounceIntegralFunction
        /** TODO: Take back if BA
        struct BounceIntegralParams {
            len_t ir; real_t xi0; real_t theta_b1; real_t theta_b2; fluxGridType fgType; 
            real_t Bmin; real_t(*F_ref)(real_t,real_t,real_t,real_t,void*); real_t(*F_eval)(real_t,real_t,real_t,real_t,void*); void *F_ref_par; int_t *Flist_eval; 
            FluxSurfaceAveragerStellarator *fsAvg; bool integrateQAWS;
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
        // this FluxSurfaceAveragerStellarator, and its generator.
        RadialGridStellarator *rGrid;
        RadialGridGeneratorStellarator *gridGenerator;

        len_t nfp = 0;

        FluxSurfaceQuantity 
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

    
        len_t nphi_interp; // number of poloidal and toroidal grid points
        real_t  
            *phi     = nullptr, // toroidal grid points
            *weights_theta = nullptr, // corresponding quadrature weights for theta
            *weights_phi   = nullptr, // corresponding quadrature weights for phi
            phi_max;

    
        void InitializeQuadrature(quadrature_method);
        void DeallocateQuadrature();
        
        static real_t FluxSurfaceIntegralFunctionThetaPhi(real_t, void *p);
        static real_t FluxSurfaceIntegralFunctionPhi(real_t, void *p);

    public:
        FluxSurfaceAveragerStellarator(
            RadialGridStellarator*, RadialGridGeneratorStellarator*, len_t nfp = 0, len_t ntheta_interp = 64, len_t nphi_interp = 64, 
            interp_method im = INTERP_LINEAR, quadrature_method qm = QUAD_FIXED_LEGENDRE
        );
        ~FluxSurfaceAveragerStellarator();

        void Rebuild();

        real_t EvaluateFluxSurfaceIntegral(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        real_t CalculateFluxSurfaceAverage(len_t ir, fluxGridType, real_t(*F)(real_t,real_t,real_t,real_t,void*), void *par=nullptr, const int_t *F_list=nullptr);
        
        const len_t GetNPhi() const
            {return nphi_interp;}    
        const real_t *GetPhi() const
            {return phi;}
        const real_t *GetWeightsTheta() const
            {return weights_theta;}
        const real_t *GetWeightsPhi() const
            {return weights_phi;}
        const real_t GetPhiMax() const    
            {return phi_max;}

        virtual FluxSurfaceQuantity *GetBOverBmin() override {return BOverBmin;}
        virtual FluxSurfaceQuantity *GetJacobian() override {return Jacobian;}
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
                return gridGenerator->BdotGradphiAtThetaPhi_f(ir,theta,phi);
            else
                return gridGenerator->BdotGradphiAtThetaPhi(ir,theta,phi);
        }
        real_t gttOverJ2AtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->gttAtThetaPhi_f(ir,theta,phi);
            else
                return gridGenerator->gttAtThetaPhi(ir,theta,phi);
        }
        real_t gtpOverJ2AtThetaPhi(len_t ir, real_t theta, real_t phi, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                return gridGenerator->gtpAtThetaPhi_f(ir,theta,phi);
            else
                return gridGenerator->gtpAtThetaPhi(ir,theta,phi);
        }
        void GeometricQuantitiesAtThetaPhi(const len_t ir, const real_t theta, const real_t phi, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2, fluxGridType fluxGridType){
            if(fluxGridType == FLUXGRIDTYPE_RADIAL)
                gridGenerator->EvaluateGeometricQuantities(ir, theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2);
            else 
                gridGenerator->EvaluateGeometricQuantities_fr(ir, theta, phi, B, Jacobian, BdotGradphi, gttOverJ2, gtpOverJ2);
        }


        len_t getNFP(){return nfp;}
        virtual const bool isStellarator() const override {return true;} 


        static real_t AssembleBAFunc(
            real_t xiOverXi0,real_t BOverBmin, real_t BdotGradphi, 
            real_t gttOverJ2, real_t gtpOverJ2, const int_t *Flist
        );
        static real_t AssembleFSAFunc(
            real_t BOverBmin, real_t BdotGradphi, real_t gttOverJ2, real_t gtpOverJ2, const int_t *Flist
        );
        
		//void PrintBOfThetaPhi(const len_t ir, const len_t N=100, enum fluxGridType fgt=FLUXGRIDTYPE_DISTRIBUTION);
    };
}

#endif/*_DREAM_FVM_FLUX_SURFACE_AVERAGER_STELLARATOR_HPP*/

