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
            const real_t *G, *G_r;            // G = R*Bphi
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
        
        real_t *r, *r_f;

        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;
        real_t diffFunc(real_t r, std::function<real_t(real_t)> F); // = dF/dr at r

        void InterpolateInputProfileToGrid(
            const len_t, const real_t*, const real_t*,
            const len_t, const real_t*,
            gsl_spline*, gsl_interp_accel*,
            real_t**, real_t**, real_t**, real_t**
        );
        gsl_spline *spline_G=nullptr, *spline_psi=nullptr, *spline_kappa=nullptr, *spline_delta=nullptr, *spline_Delta=nullptr;
        gsl_interp_accel *gsl_acc_G, *gsl_acc_psi, *gsl_acc_kappa, *gsl_acc_delta, *gsl_acc_Delta;

        real_t normalizedJacobian(len_t,real_t);
        real_t normalizedJacobian(len_t,real_t,real_t,real_t);
        real_t normalizedJacobian_f(len_t,real_t);
        real_t normalizedJacobian_f(len_t,real_t,real_t,real_t);
    
    public:
        AnalyticBRadialGridGenerator(
            const len_t nr, real_t r0, real_t ra, real_t R0,
            len_t ntheta_interp, struct shape_profiles*
        );
        ~AnalyticBRadialGridGenerator();

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void DeallocateShapeProfiles();

        virtual real_t JacobianAtTheta(const len_t ir, const real_t) override;
        virtual real_t JacobianAtTheta(const len_t ir, const real_t, const real_t, const real_t) override;
        virtual real_t ROverR0AtTheta(const len_t, const real_t) override;
        virtual real_t ROverR0AtTheta(const len_t, const real_t, const real_t, const real_t) override;
        virtual real_t NablaR2AtTheta(const len_t, const real_t) override;
        virtual real_t NablaR2AtTheta(const len_t, const real_t, const real_t, const real_t) override;
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t) override;
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t, const real_t, const real_t) override;
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t) override;
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t, const real_t, const real_t) override;
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t) override;
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t, const real_t, const real_t) override;
        virtual void EvaluateGeometricQuantities(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override;
        virtual void EvaluateGeometricQuantities_fr(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override;

    };
}

#endif/*_DREAM_FVM_GRID_ANALYTIC_B_RADIAL_GRID_GENERATOR_HPP*/
