#ifndef _DREAM_EMPTY_RADIAL_GRID_HPP
#define _DREAM_EMPTY_RADIAL_GRID_HPP

#include "FVM/Grid/RadialGridGenerator.hpp"

namespace DREAM::FVM {
    class EmptyRadialGridGenerator : public RadialGridGenerator {
        public:
            EmptyRadialGridGenerator() : RadialGridGenerator(1) {
                ntheta_interp = 1;
                isUpDownSymmetric = true;
            }
            bool isBuilt = false;

            virtual bool NeedsRebuild(const real_t) const override { return !isBuilt; }
            virtual bool Rebuild(const real_t, RadialGrid*) override;
            virtual real_t JacobianAtTheta(const len_t, const real_t) override {return 1.0;}
            virtual real_t ROverR0AtTheta(const len_t, const real_t) override {return 1.0;}
            virtual real_t NablaR2AtTheta(const len_t, const real_t) override {return 1.0;}
            virtual real_t JacobianAtTheta_f(const len_t, const real_t) override {return 1.0;}
            virtual real_t ROverR0AtTheta_f(const len_t, const real_t) override {return 1.0;}
            virtual real_t NablaR2AtTheta_f(const len_t, const real_t) override {return 1.0;}
            virtual void EvaluateGeometricQuantities(const len_t, const real_t, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override
                {Jacobian=1; B=1; NablaR2 = 1; ROverR0 = 1;}
            virtual void EvaluateGeometricQuantities_fr(const len_t, const real_t, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) override
                {Jacobian=1; B=1; NablaR2 = 1; ROverR0 = 1;}

			virtual void GetRThetaPhiFromCartesian(real_t*, real_t*, real_t*, real_t, real_t , real_t , real_t, real_t ) override {}
			virtual void GetGradRCartesian(real_t*, real_t , real_t, real_t) override {}

			virtual const real_t GetZ0() override { return 0; }
			virtual const len_t GetNPsi() override { return 0; }
			virtual const len_t GetNTheta() override { return 0; }
			virtual const real_t *GetFluxSurfaceRMinusR0() override { return nullptr; }
			virtual const real_t *GetFluxSurfaceRMinusR0_f() override { return nullptr; }
			virtual const real_t *GetFluxSurfaceZMinusZ0() override { return nullptr; }
			virtual const real_t *GetFluxSurfaceZMinusZ0_f() override { return nullptr; }
			virtual const real_t *GetPoloidalAngle() override { return nullptr; }
    };

    class EmptyRadialGrid : public RadialGrid {
        public:
            EmptyRadialGrid(): RadialGrid(new EmptyRadialGridGenerator()) {}

            virtual ~EmptyRadialGrid() {}

    };
}

#endif/*_DREAM_EMPTY_RADIAL_GRID_HPP*/
