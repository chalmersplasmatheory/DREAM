#ifndef _DREAM_FVM_MOMENTUM_GRID_GENERATOR_HPP
#define _DREAM_FVM_MOMENTUM_GRID_GENERATOR_HPP

#include <string>
#include <utility>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

namespace DREAM::FVM {
    class MomentumGridGenerator {
    public:
        virtual ~MomentumGridGenerator() {}

        virtual bool NeedsRebuild(const real_t, const bool) = 0;
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*) = 0;
    };

    class MomentumGridGeneratorException : public FVMException {
    public:
        template<typename ... Args>
        MomentumGridGeneratorException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("MomentumGridGenerator");
        }
    };
}

#endif/*_DREAM_FVM_MOMENTUM_GRID_GENERATOR_HPP*/
