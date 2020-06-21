//namespace DREAM::FVM { class RadialGridGenerator; }

#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include <functional>
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"

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
        len_t ntheta_interp;
        len_t ntheta_ref;
        real_t *theta_ref;
        real_t R0;
        real_t **B_ref = nullptr, **Jacobian_ref,
               **ROverR0_ref,     **NablaR2_ref,
               **B_ref_f,         **Jacobian_ref_f,
               **ROverR0_ref_f,   **NablaR2_ref_f,
                *Bmin,             *Bmin_f,
                *Bmax,             *Bmax_f,
                *BtorGOverR0,      *BtorGOverR0_f;

        // True if the flux surfaces are up-down symmetric, i.e. if B(theta) = B(-theta)
        // where (if true) theta=0 must correspond to outermost low-field side, B(0) = B_min. 
        bool isUpDownSymmetric = false;

        void SetNr(const len_t n) { this->nr = n; }

    public:

        RadialGridGenerator(const len_t nr); 
        virtual ~RadialGridGenerator();
        
        len_t GetNr() const { return this->nr; }
        len_t GetNthetaInterp() const { return this->ntheta_interp; }

        virtual bool NeedsRebuild(const real_t t) const = 0;
        virtual bool Rebuild(const real_t t, RadialGrid*) = 0;
        virtual void RebuildJacobians(RadialGrid*);
        virtual void CreateMagneticFieldData(const real_t *x, const real_t *x_f) = 0;
        bool IsFieldSymmetric(){return isUpDownSymmetric;}
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_GENERATOR_HPP*/
