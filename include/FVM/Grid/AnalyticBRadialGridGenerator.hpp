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
    public:
        struct shape_profiles {
            len_t nG, npsi, nkappa, ndelta, nDelta;
            const real_t *GOverR0, *G_r;      // G/R0 = R/R0*Bphi
            const real_t *psi, *psi_r;        // Poloidal flux
            const real_t *kappa, *kappa_r;    // Elongation
            const real_t *delta, *delta_r;    // Triangularity
            const real_t *Delta, *Delta_r;    // Shafranov shift
        };

    private:
        real_t rMin, rMax;
        struct shape_profiles *providedProfiles;

        real_t *psi = nullptr, *kappa, *delta, *Delta,
            *GPrime, *kappaPrime, *deltaPrime, *DeltaPrime;
        real_t *psi_f, *kappa_f, *delta_f, *Delta_f,
            *GPrime_f, *kappaPrime_f, *deltaPrime_f, *DeltaPrime_f;
        
        real_t *rf_provided=nullptr;
        real_t *r, *r_f=nullptr;

        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;

        bool R0IsInf;

        void InterpolateInputProfileToGrid(
            const len_t, const real_t*, const real_t*,
            const len_t, const real_t*, const real_t*,
            gsl_spline*, gsl_interp_accel*,
            real_t**, real_t**, real_t**, real_t**
        );
        
        real_t InterpolateInputProfileSingleExtrap(real_t r,
    		const len_t nProvided, const real_t *xProvided, const real_t *xProvided_r,
    		gsl_spline *spline_x, gsl_interp_accel *spline_acc
			);
        real_t InterpolateInputProfileSingleDerivExtrap(real_t r,
    		const len_t nProvided, const real_t *xProvided_r,
    		gsl_spline *spline_x, gsl_interp_accel *spline_acc
			);
		real_t InterpolateInputElongation(real_t r);
		real_t InterpolateInputElongationDeriv(real_t r);
		real_t InterpolateInputTriangularity(real_t r);
		real_t InterpolateInputTriangularityDeriv(real_t r);
		real_t InterpolateInputShafranovShift(real_t r);
		real_t InterpolateInputShafranovShiftDeriv(real_t r);
			
        gsl_spline *spline_G=nullptr, *spline_psi=nullptr, *spline_kappa=nullptr, *spline_delta=nullptr, *spline_Delta=nullptr;
        gsl_interp_accel *gsl_acc_G, *gsl_acc_psi, *gsl_acc_kappa, *gsl_acc_delta, *gsl_acc_Delta;

        real_t normalizedJacobian(len_t ir,real_t theta) 
            {return normalizedJacobian(ir,theta,cos(theta),sin(theta));}
        real_t normalizedJacobian(len_t,real_t,real_t,real_t);
        real_t normalizedJacobian_f(len_t ir,real_t theta)
            {return normalizedJacobian_f(ir,theta,cos(theta),sin(theta));}
        real_t normalizedJacobian_f(len_t,real_t,real_t,real_t);
    
        void constructSplines(struct shape_profiles*);
    public:
        AnalyticBRadialGridGenerator(
            const len_t nr, real_t r0, real_t ra, real_t R0,
            len_t ntheta_interp, struct shape_profiles*
        );
        AnalyticBRadialGridGenerator(
            const real_t *r_f, const len_t nr, real_t R0,
            len_t ntheta_interp, struct shape_profiles*
        );
        ~AnalyticBRadialGridGenerator();

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void DeallocateShapeProfiles();

        virtual real_t JacobianAtTheta(const len_t ir, const real_t theta) override;
        virtual real_t ROverR0AtTheta(const len_t ir, const real_t theta) override;
        virtual real_t NablaR2AtTheta(const len_t ir, const real_t theta) override;
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t theta) override;
        virtual real_t ROverR0AtTheta_f(const len_t ir, const real_t theta) override;
        virtual real_t NablaR2AtTheta_f(const len_t ir, const real_t theta) override;
        virtual void EvaluateGeometricQuantities(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override;
        virtual void EvaluateGeometricQuantities_fr(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override;
        
		virtual void GetRThetaPhiFromCartesian(real_t*, real_t*, real_t*, real_t, real_t, real_t, real_t, real_t ) override;
		virtual void GetGradRCartesian(real_t*, real_t , real_t, real_t ) override;

		virtual const real_t GetZ0() override { return 0; }
		virtual const len_t GetNPsi() override { return this->GetNr(); }
		virtual const len_t GetNTheta() override { return 120; }
		virtual const real_t *GetFluxSurfaceROverR0() override;
		virtual const real_t *GetFluxSurfaceROverR0_f() override;
		virtual const real_t *GetFluxSurfaceZ() override;
		virtual const real_t *GetFluxSurfaceZ_f() override;
		virtual const real_t *GetPoloidalAngle() override;
    };
}

#endif/*_DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP*/
