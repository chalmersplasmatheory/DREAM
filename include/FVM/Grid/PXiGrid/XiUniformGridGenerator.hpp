#ifndef _TQS_FVM_P_XI_GRID_XI_UNIFORM_GRID_GENERATOR_HPP
#define _TQS_FVM_P_XI_GRID_XI_UNIFORM_GRID_GENERATOR_HPP

#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Grid/PXiGrid/XiGridGenerator.hpp"

namespace TQS::FVM::PXiGrid {
    class XiUniformGridGenerator : public XiGridGenerator {
    private:
        len_t nxi;
        real_t xiMin, xiMax;

        bool initialized = false;
    public:
        XiUniformGridGenerator(const len_t);

        virtual bool NeedsRebuild(const real_t, const bool) { return (!initialized); }
        virtual bool Rebuild(const real_t, const len_t, MomentumGrid*, const RadialGrid*);
    };
}

#endif/*_TQS_FVM_P_XI_GRID_XI_UNIFORM_GRID_GENERATOR_HPP*/
