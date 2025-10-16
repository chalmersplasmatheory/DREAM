//namespace DREAM::FVM { class RadialGridGeneratorStellarator; }

#include "FVM/Grid/Stellarator/RadialGridStellarator.hpp"
#include "gsl/gsl_multimin.h"

#ifndef _DREAM_FVM_RADIAL_GRID_GENERATOR_STELLARATOR_HPP
#define _DREAM_FVM_RADIAL_GRID_GENERATOR_STELLARATOR_HPP

#include "FVM/Grid/RadialGridGenerator.hpp"

/***************************************************
 * Base class for radial grid generators. *
 ***************************************************/

namespace DREAM::FVM {
    class RadialGridGeneratorStellarator : public RadialGridGenerator {

    private:
        len_t nr = 0;

    protected:
        virtual real_t FindMagneticFieldExtremum(
			len_t ir, int_t sgn, enum fluxGridType
		) override;
        gsl_multimin_fminimizer *gsl_multi_fmin;

        len_t nphi_interp;
        real_t  *BpolIOverR0,        *BpolIOverR0_f,        // poloidal magnetic field strength (I(r)/R0)
                *iota,               *iota_f,               // Rotational transform
                *Fpassing,           *Fpassing_f;
       
        // Tolerance to use when switching between cartesian and flux surface coordinates
        // Increasing the value of this tolerance seem to make it hard for SPI simulations
        // to converge without also increasing the overall generak tolerance for the Newton solver,
        // which is why we use this rather small value here. It does however not seem to slow down the
        // simulations significantly
        
        len_t nfp;  // Number of field poles?
        
    public:
        RadialGridGeneratorStellarator(const len_t nr); 
        virtual ~RadialGridGeneratorStellarator();
        
        len_t GetNphiInterp() const { return this->nphi_interp; }

        virtual bool NeedsRebuild(const real_t t) const = 0;
        virtual bool Rebuild(const real_t , RadialGrid*) {return false;}
        virtual bool Rebuild(const real_t t, RadialGridStellarator*) = 0;
        virtual void RebuildJacobians(RadialGridStellarator*);
        len_t getNFP(){return nfp;}

        virtual real_t BAtThetaPhi(const len_t, const real_t, const real_t) = 0;
        virtual real_t BAtThetaPhi_f(const len_t, const real_t, const real_t) = 0;
        
        // The following functions set the geometry and are implemented in derived classes
        virtual real_t JacobianAtThetaPhi(const len_t, const real_t, const real_t) = 0;
        virtual real_t JacobianAtThetaPhi_f(const len_t, const real_t, const real_t) = 0;
        virtual real_t BdotGradphiAtThetaPhi(const len_t, const real_t, const real_t) = 0;
        virtual real_t BdotGradphiAtThetaPhi_f(const len_t, const real_t, const real_t) = 0;
        virtual real_t gttAtThetaPhi(const len_t, const real_t, const real_t) = 0;
        virtual real_t gttAtThetaPhi_f(const len_t, const real_t, const real_t) = 0;
        virtual real_t gtpAtThetaPhi(const len_t, const real_t, const real_t) = 0;
        virtual real_t gtpAtThetaPhi_f(const len_t, const real_t, const real_t) = 0;
        
        virtual void EvaluateGeometricQuantities(const len_t, const real_t, const real_t, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2) = 0;
        virtual void EvaluateGeometricQuantities_fr(const len_t, const real_t, const real_t, real_t &B, real_t &Jacobian, real_t &BdotGradphi, real_t &gttOverJ2, real_t &gtpOverJ2) = 0;

		// Helper routines for saving equilibrium to output
		virtual const len_t GetNPhi() = 0;
		virtual const real_t *GetPoloidalAngle() = 0;
		virtual const real_t *GetToroidalAngle() = 0;


        // Functions that are unused for a stellarator 
        // TODO: Is this ok?

        virtual void GetRThetaPhiFromCartesian(real_t*, real_t*, real_t*, real_t , real_t , real_t, real_t, real_t ){};
        virtual void GetGradRCartesian(real_t* ,real_t, real_t, real_t){};
        virtual real_t FindClosestApproach(real_t , real_t , real_t , real_t , real_t , real_t ){return 0;};
        
        virtual real_t BAtTheta(const len_t, const real_t) {return 0;}
        virtual real_t BAtTheta_f(const len_t, const real_t) {return 0;}
        
        virtual real_t JacobianAtTheta(const len_t, const real_t) {return 0;}
        virtual real_t ROverR0AtTheta(const len_t, const real_t) {return 0;}
        virtual real_t NablaR2AtTheta(const len_t, const real_t) {return 0;}
        virtual real_t JacobianAtTheta_f(const len_t, const real_t) {return 0;}
        virtual real_t ROverR0AtTheta_f(const len_t, const real_t) {return 0;}
        virtual real_t NablaR2AtTheta_f(const len_t, const real_t) {return 0;}

        virtual void EvaluateGeometricQuantities(const len_t, const real_t, real_t &, real_t &, real_t &, real_t &) {}
        virtual void EvaluateGeometricQuantities_fr(const len_t, const real_t, real_t &, real_t &, real_t &, real_t &) {}

        virtual const real_t GetZ0() {return 0;}
		virtual const real_t *GetFluxSurfaceRMinusR0() {return nullptr;}
		virtual const real_t *GetFluxSurfaceRMinusR0_f() {return nullptr;}
		virtual const real_t *GetFluxSurfaceZMinusZ0() {return nullptr;}
		virtual const real_t *GetFluxSurfaceZMinusZ0_f() {return nullptr;}
        
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_GENERATOR_STELLARATOR_HPP*/
