#ifndef _DREAMTESTS_KNOCKON_TEST_GRID_HPP
#define _DREAMTESTS_KNOCKON_TEST_GRID_HPP

#include <memory>
#include <vector>

#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiCustomGridGenerator.hpp"

namespace DREAMTESTS::_DREAM::KnockOnTest {
// Grid owns its momentum grids, but the radial grid has a separate lifetime.
struct GridDeleter {
    void operator()(DREAM::FVM::Grid *grid) const {
        auto *radial = grid->GetRadialGrid();
        delete grid;
        delete radial;
    }
};
using OwnedGrid = std::unique_ptr<DREAM::FVM::Grid, GridDeleter>;

inline std::vector<real_t> UniformPitchFaces(len_t n) {
    std::vector<real_t> faces(n+1);
    for (len_t j=0; j<=n; ++j) faces[j] = -1 + 2.0*j/n;
    return faces;
}

// A narrow annulus about r/R0=epsilon. Constant profile arrays outlive the
// generator, which owns the heap-allocated descriptor. No production-sized
// momentum grid is needed to prepare a genuine native toroidal geometry.
inline OwnedGrid ToroidalGrid(
    real_t epsilon, const std::vector<real_t>& faces, len_t nr=1,
    len_t np=1, real_t pMin=0, real_t pMax=2
) {
    using namespace DREAM::FVM;
    static const real_t r[] = {0,1}, one[] = {1,1}, zero[] = {0,0};
    auto *profiles = new AnalyticBRadialGridGenerator::shape_profiles{
        2,2,2,2,2, one,r,zero,r,one,r,zero,r,zero,r
    };
    auto *generator = new AnalyticBRadialGridGenerator(
        nr,epsilon-0.01,epsilon+0.01,1,128,profiles
    );
    auto *radial = new RadialGrid(generator);
    auto *momentum = new PXiGrid::PXiMomentumGrid(new PXiGrid::MomentumGridGenerator(
        new PXiGrid::PUniformGridGenerator(np,pMin,pMax),
        new PXiGrid::XiCustomGridGenerator(faces.data(),faces.size()-1)
    ),0,radial);
    OwnedGrid grid(new Grid(radial,momentum));
    grid->Rebuild(0);
    return grid;
}
} // namespace DREAMTESTS::_DREAM::KnockOnTest
#endif
