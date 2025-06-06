//namespace DREAM::FVM { class RadialGridGenerator; }

#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include <functional>
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_min.h"

#ifndef _DREAM_FVM_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_RADIAL_GRID_GENERATOR_HPP

/***************************************************
 * Base class for radial grid generators. *
 ***************************************************/

namespace DREAM::FVM {
    class RadialGridGenerator {

    private:
        len_t nr = 0;

    protected:
        virtual real_t FindMagneticFieldExtremum(
			len_t ir, int_t sgn, enum fluxGridType
		);
        gsl_min_fminimizer *gsl_fmin;

        len_t ntheta_interp;
        real_t R0;
        real_t  *Bmin=nullptr,       *Bmin_f,               // minimum B on flux surface
                *Bmax,               *Bmax_f,               // maximum B on flux surface
                *theta_Bmin,         *theta_Bmin_f,         // poloidal angle of minimum B
                *theta_Bmax,         *theta_Bmax_f,         // poloidal angle of maximum B
                *xi0TrappedBoundary, *xi0TrappedBoundary_f, // pitch of trapped-passing boundary
                *BtorGOverR0,        *BtorGOverR0_f,        // toroidal magnetic field strength (G(r)/R0, where G(r)/R is the true toroidal field) 
                *psiPrimeRef,        *psiPrimeRef_f;        // radial derivative of reference poloidal flux (normalized to R0), (in cylindrical geometry is the poloidal magnetic field strength)
        // TODO: in the interface, specify jtotRef instead of psiPrimeRef, where jtotRef
        //       is a reference current profile from which we calculate the generated poloidal magnetic field (or flux derivative) 

        // Tolerance to use when switching between cartesian and flux surface coordinates
        // Increasing the value of this tolerance seem to make it hard for SPI simulations
        // to converge without also increasing the overall generak tolerance for the Newton solver,
        // which is why we use this rather small value here. It does however not seem to slow down the
        // simulations significantly
        const real_t CartesianCoordinateTol = 1e-6;
        
        // True if the flux surfaces are up-down symmetric, i.e. if B(theta) = B(-theta)
        // where (if true) theta=0 must correspond to outermost low-field side, B(0) = B_min. 
        bool isUpDownSymmetric = false;

        void SetNr(const len_t n) { this->nr = n; }

        // Return poloidal angle of minimum and maximum magnetic fields.
        // In principle, for extreme poloidal fields, the global minimum will not occur at 0,
        // and we'll handle that scenario with the below functions.
        virtual real_t getTheta_Bmin(const len_t ir) {
            return FindMagneticFieldExtremum(ir,1,FLUXGRIDTYPE_DISTRIBUTION);
        }
        virtual real_t getTheta_Bmax(const len_t ir) {
            return FindMagneticFieldExtremum(ir,-1,FLUXGRIDTYPE_DISTRIBUTION);
        }
        virtual real_t getTheta_Bmin_f(const len_t ir) {
            return FindMagneticFieldExtremum(ir,1,FLUXGRIDTYPE_RADIAL);
        }
        virtual real_t getTheta_Bmax_f(const len_t ir) {
            return FindMagneticFieldExtremum(ir,-1,FLUXGRIDTYPE_RADIAL);
        }
        
    public:
        RadialGridGenerator(const len_t nr); 
        virtual ~RadialGridGenerator();
        
        len_t GetNr() const { return this->nr; }
        len_t GetNthetaInterp() const { return this->ntheta_interp; }

        virtual void GetRThetaPhiFromCartesian(real_t*, real_t*, real_t*, real_t , real_t , real_t, real_t, real_t )=0;
        virtual void GetGradRCartesian(real_t* ,real_t, real_t, real_t)=0;
        virtual real_t FindClosestApproach(real_t , real_t , real_t , real_t , real_t , real_t );

        virtual bool NeedsRebuild(const real_t t) const = 0;
        virtual bool Rebuild(const real_t t, RadialGrid*) = 0;
        virtual void RebuildJacobians(RadialGrid*);
        bool IsFieldSymmetric(){return isUpDownSymmetric;}

        virtual real_t BAtTheta(const len_t ir, const real_t theta);
        virtual real_t BAtTheta_f(const len_t ir, const real_t theta);
        
        // The following functions set the geometry and are implemented in derived classes
        virtual real_t JacobianAtTheta(const len_t ir, const real_t theta) = 0;
        virtual real_t ROverR0AtTheta(const len_t ir, const real_t theta) = 0;
        virtual real_t NablaR2AtTheta(const len_t ir, const real_t theta) = 0;
        virtual real_t JacobianAtTheta_f(const len_t ir, const real_t theta) = 0;
        virtual real_t ROverR0AtTheta_f(const len_t ir, const real_t theta) = 0;
        virtual real_t NablaR2AtTheta_f(const len_t ir, const real_t theta) = 0;
        
        virtual void EvaluateGeometricQuantities(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) = 0;
        virtual void EvaluateGeometricQuantities_fr(const len_t ir, const real_t theta, real_t &B, real_t &Jacobian, real_t &ROverR0, real_t &NablaR2) = 0;

		// Helper routines for saving equilibrium to output
		virtual const real_t GetZ0() = 0;
		virtual const len_t GetNPsi() = 0;
		virtual const len_t GetNTheta() = 0;
		virtual const real_t *GetFluxSurfaceRMinusR0() = 0;
		virtual const real_t *GetFluxSurfaceRMinusR0_f() = 0;
		virtual const real_t *GetFluxSurfaceZMinusZ0() = 0;
		virtual const real_t *GetFluxSurfaceZMinusZ0_f() = 0;
		virtual const real_t *GetPoloidalAngle() = 0;
        
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_GENERATOR_HPP*/
