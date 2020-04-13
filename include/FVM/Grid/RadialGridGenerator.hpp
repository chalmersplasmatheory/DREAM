namespace DREAM::FVM { class RadialGridGenerator; }

#include "FVM/config.h"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include <functional>
#include "FVM/Grid/MagneticQuantityHandler.hpp"

#ifndef _DREAM_FVM_RADIAL_GRID_GENERATOR_HPP
#define _DREAM_FVM_RADIAL_GRID_GENERATOR_HPP

/***************************************************
 * Abstract base class for radial grid generators. *
 ***************************************************/

namespace DREAM::FVM {
    class RadialGridGenerator {
    private:
        len_t nr=0;

    protected:
        void SetNr(const len_t n) { this->nr = n; }

    public:
        RadialGridGenerator(const len_t nr) : nr(nr) {}
        virtual ~RadialGridGenerator() {};
        
        len_t GetNr() const { return this->nr; }

        virtual bool NeedsRebuild(const real_t t) const = 0;
        virtual bool Rebuild(const real_t t, RadialGrid*) = 0;
        virtual void RebuildJacobians(RadialGrid*, MomentumGrid**,MagneticQuantityHandler*) = 0;
    };
}

#endif/*_DREAM_FVM_RADIAL_GRID_GENERATOR_HPP*/
