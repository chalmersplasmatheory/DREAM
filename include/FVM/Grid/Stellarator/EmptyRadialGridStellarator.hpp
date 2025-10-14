#ifndef _DREAM_EMPTY_RADIAL_GRID_HPP // TODO
#define _DREAM_EMPTY_RADIAL_GRID_HPP

#include "FVM/Grid/RadialGridGenerator.hpp"

namespace DREAM::FVM {
    class EmptyRadialGridGenerator : public RadialGridGenerator {
        protected:
            len_t nphi_interp;
        public:
            EmptyRadialGridGenerator() : RadialGridGenerator(1) {
                ntheta_interp = 1;
                nphi_interp = 1;
                isUpDownSymmetric = true;
            }
            bool isBuilt = false;

            virtual bool NeedsRebuild(const real_t) const override { return !isBuilt; }
            virtual bool Rebuild(const real_t, RadialGrid*) override;
            virtual real_t JacobianAtThetaPhi(const len_t, const real_t, const real_t) override {return 1.0;}
            virtual real_t BdotGradphiAtThetaPhi(const len_t, const real_t, const real_t)  {return 1.0;}
            virtual real_t gttAtThetaPhi(const len_t, const real_t, const real_t)  {return 1.0;}
            virtual real_t gtpAtThetaPhi(const len_t, const real_t, const real_t)  {return 1.0;}
            virtual real_t JacobianAtThetaPhi_f(const len_t, const real_t, const real_t) override {return 1.0;}
            virtual real_t BdotGradphiAtThetaPhi_f(const len_t, const real_t, const real_t)  {return 1.0;}
            virtual real_t gttAtThetaPhi_f(const len_t, const real_t, const real_t)  {return 1.0;}
            virtual real_t gtpAtThetaPhi_f(const len_t, const real_t, const real_t)  {return 1.0;}
            
            virtual void EvaluateGeometricQuantities(const len_t, const real_t, const real_t, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2) 
                {B=1; Jacobian=1; BdotGradphi=1; gttOverJ2=1; gtpOverJ2=1}
            virtual void EvaluateGeometricQuantities_fr(const len_t, const real_t, const real_t, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2) 
                {B=1; Jacobian=1; BdotGradphi=1; gttOverJ2=1; gtpOverJ2=1}

			virtual const len_t GetNPsi() override { return 0; }
			virtual const len_t GetNTheta() override { return 0; }
            virtual const len_t GetNPhi() override { return 0; }
			virtual const real_t *GetPoloidalAngle() override { return nullptr; }
            virtual const real_t *GetToroidalAngle() override { return nullptr; }
    };

    class EmptyRadialGrid : public RadialGrid {
        public:
            EmptyRadialGrid(): RadialGrid(new EmptyRadialGridGenerator()) {}

            virtual ~EmptyRadialGrid() {}

    };
}

#endif/*_DREAM_EMPTY_RADIAL_GRID_HPP*/
