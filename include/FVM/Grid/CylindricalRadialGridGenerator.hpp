#ifndef _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/RadialGridGenerator.hpp"
#include "FVM/Grid/MagneticQuantityHandler.hpp"
#include <functional>

namespace DREAM::FVM {
    class CylindricalRadialGridGenerator : public RadialGridGenerator {
    private:
        real_t xMin=0, xMax=1;
        real_t B0=0;
        
        len_t ntheta_ref;
        real_t *theta_ref;
        real_t **B_ref,         **Jacobian_ref,
               **ROverR0_ref,   **NablaR2_ref,
               **B_ref_f,       **Jacobian_ref_f,
               **ROverR0_ref_f, **NablaR2_ref_f,
                *Bmin_ref,       *Bmin_ref_f,
                *Bmax_ref,       *Bmax_ref_f;

        // Set to true when the grid is constructed for the first time
        bool isBuilt = false;

    public:
        CylindricalRadialGridGenerator(const len_t nx, const real_t B0, const real_t x0=0, const real_t xa=1);

        virtual bool NeedsRebuild(const real_t) const override { return (!isBuilt); }
        virtual bool Rebuild(const real_t, RadialGrid*) override;
        virtual void CreateMagneticFieldData(const real_t *x, const real_t *x_f);
        virtual void RebuildJacobians(RadialGrid*, MomentumGrid**,MagneticQuantityHandler*) override;
        virtual void DeallocateMagneticFieldData();

    };
}

#endif/*_DREAM_FVM_GRID_CYLINDRICAL_RADIAL_GRID_GENERATOR_HPP*/
