#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"

using namespace DREAM;

MollerDeltaAngleKernel::MollerDeltaAngleKernel(
    const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t p_cutoff,
    len_t n_xi_stars_tabulate, len_t poloidal_order, len_t pitch_order
)
    : gridK(grid_knockon),
      gridP(grid_primary),
      pCutoff(p_cutoff),
      nXiStarsTabulate(n_xi_stars_tabulate),
      poloidalOrder(poloidal_order), pitchOrder(pitch_order) {
    ValidateInputParameters();
    GridRebuilt();
}

MollerDeltaAngleKernel::~MollerDeltaAngleKernel() { Deallocate(); }

void MollerDeltaAngleKernel::Deallocate() {
    const len_t Nr = nrAllocated;

    if (deltaTable != nullptr) {
        for (len_t ir = 0; ir < Nr; ++ir) {
            delete[] deltaTableStorage[ir];
            delete[] deltaTable[ir];
        }
        delete[] deltaTableStorage;
        delete[] deltaTable;
        deltaTableStorage = nullptr;
        deltaTable = nullptr;
    }

    if (xiInterp != nullptr) {
        for (len_t ir = 0; ir < Nr; ++ir) delete[] xiInterp[ir];
        delete[] xiInterp;
        xiInterp = nullptr;
    }

    if (xiStarsTab != nullptr) {
        delete[] xiStarsTab;
        xiStarsTab = nullptr;
    }
    nrAllocated = 0;
}

void DREAM::MollerDeltaAngleKernel::ValidateGridAssumptions() const {
    if (gridK == nullptr || gridP == nullptr)
        throw DREAMException("MollerDeltaAngleKernel: grid pointers must not be null.");

    if (!std::isfinite(pCutoff) || !(pCutoff > 0))
        throw DREAMException(
            "MollerDeltaAngleKernel: invalid pCutoff=%.16g (must be > 0).", pCutoff
        );

    if (nXiStarsTabulate < 2)
        throw DREAMException(
            "MollerDeltaAngleKernel: nXiStarsTabulate must be >= 2 (got %ld).",
            (long)nXiStarsTabulate
        );
    if (poloidalOrder < 1)
        throw DREAMException(
            "MollerDeltaAngleKernel: poloidalOrder must be >= 1 (got %ld).", (long)poloidalOrder
        );

    const len_t NrK = gridK->GetNr();
    const len_t NrP = gridP->GetNr();
    if (NrK != NrP)
        throw DREAMException("MollerDeltaAngleKernel: gridK and gridP must have same Nr.");

    for (len_t ir = 1; ir < NrK; ++ir) {
        if (gridK->GetNp1(ir) != gridK->GetNp1(0) || gridK->GetNp2(ir) != gridK->GetNp2(0)) {
            throw DREAMException(
                "MollerDeltaAngleKernel: requires uniform knock-on momentum grid across radii."
            );
        }
        if (gridP->GetNp1(ir) != gridP->GetNp1(0) || gridP->GetNp2(ir) != gridP->GetNp2(0)) {
            throw DREAMException(
                "MollerDeltaAngleKernel: requires uniform primary momentum grid across radii."
            );
        }
    }
}

void DREAM::MollerDeltaAngleKernel::ValidateInputParameters() const {
    if (pitchOrder<1)
        throw DREAMException("MollerDeltaAngleKernel: pitchOrder must be >= 1.");

    if (nXiStarsTabulate < 2) {
        throw DREAMException(
            "MollerDeltaAngleKernel: nXiStarsTabulate must be >= 2 (got %ld).",
            (long)nXiStarsTabulate
        );
    }
    if (poloidalOrder < 1) {
        throw DREAMException(
            "MollerDeltaAngleKernel: poloidalOrder must be >= 1 (got %ld).", (long)poloidalOrder
        );
    }
}

// Use exactly the same non-zero momentum-cell support as the energy kernel.
// In a partially accessible cell, keep the representative momentum inside
// the integration interval: a coarse cell centre can lie above p1 and yield
// xi_star > 1 even though part of the cell has a non-zero cross section.
bool MollerDeltaAngleKernel::EvaluateXiStarOnGrids(
    len_t i, len_t k, real_t &xiStar
) const {
    if (KnockOnUtilities::EvaluateMollerFluxMatrixElementOnGrid(
            i, k, gridK, gridP, pCutoff) == 0)
        return false;

    const auto *mgK = gridK->GetMomentumGrid(0);
    const auto *mgP = gridP->GetMomentumGrid(0);
    const real_t p1 = mgP->GetP1(k);
    const real_t pLower = std::max(pCutoff, mgK->GetP1_f(i));
    const real_t pUpper = std::min(
        mgK->GetP1_f(i+1), KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1)
    );
    const real_t p = std::max(pLower, std::min(mgK->GetP1(i), pUpper));
    xiStar = KnockOnUtilities::Kinematics::EvaluateXiStar(p, p1);
    if (!std::isfinite(xiStar) || xiStar < 0 || xiStar > 1)
        throw DREAMException(
            "MollerDeltaAngleKernel: invalid xi_star=%.16g at i=%ld, k=%ld.",
            xiStar, (long)i, (long)k
        );
    return true;
}

// Evaluate the range of angles which will actually be used by the operator.
void MollerDeltaAngleKernel::EvaluateXiStarRangeOnGrids(real_t &xiStarMin, real_t &xiStarMax) {
    auto *mgK = gridK->GetMomentumGrid(0);
    auto *mgP = gridP->GetMomentumGrid(0);
    xiStarMin = +std::numeric_limits<real_t>::infinity();
    xiStarMax = -std::numeric_limits<real_t>::infinity();
    for (len_t i = 0; i < mgK->GetNp1(); ++i) {
        for (len_t k = 0; k < mgP->GetNp1(); ++k) {
            real_t xs;
            if (!EvaluateXiStarOnGrids(i, k, xs)) {
                continue;
            }
            xiStarMin = std::min(xiStarMin, xs);
            xiStarMax = std::max(xiStarMax, xs);
        }
    }
}

// Tabulate the kinematic delta over the range of xi_star reached by
// momentum-cell pairs that contribute to the knock-on source.
void MollerDeltaAngleKernel::BuildXiStarTabulationGrid() {
    EvaluateXiStarRangeOnGrids(xiStarMin, xiStarMax);
    hasActivePairs = std::isfinite(xiStarMin) && std::isfinite(xiStarMax);
    // A grid with no allowed collisions is a valid zero source.
    if (!hasActivePairs)
        xiStarMin = xiStarMax = 0;
    dXiStar = (xiStarMax - xiStarMin) / (real_t)(nXiStarsTabulate - 1);

    xiStarsTab = new real_t[nXiStarsTabulate];
    for (len_t i = 0; i < nXiStarsTabulate; ++i) {
        xiStarsTab[i] = xiStarMin + (real_t)i * dXiStar;
    }
}

void MollerDeltaAngleKernel::AllocateDeltaTables() {
    const len_t Nr = gridK->GetNr();
    nrAllocated = Nr;

    deltaTable = new real_t **[Nr];
    deltaTableStorage = new real_t *[Nr];

    for (len_t ir = 0; ir < Nr; ++ir) {
        const len_t NxiK = gridK->GetNp2(ir);
        const len_t NxiP = gridP->GetNp2(ir);
        const len_t sliceSize = NxiK * NxiP;

        deltaTable[ir] = new real_t *[nXiStarsTabulate];
        deltaTableStorage[ir] = new real_t[nXiStarsTabulate * sliceSize]();

        for (len_t m = 0; m < nXiStarsTabulate; ++m) {
            deltaTable[ir][m] = deltaTableStorage[ir] + m * sliceSize;
        }
    }

    xiInterp = new XiStarInterp *[Nr];
    for (len_t ir = 0; ir < Nr; ++ir) {
        xiInterp[ir] = new XiStarInterp[gridK->GetNp1(ir) * gridP->GetNp1(ir)];
    }
}

/**
 * Tabulate the pitch redistribution as a function of the scattering cosine
 * xi_star(p, p1), which carries the momentum dependence of the angular kernel.
 * Each column integrates one primary pitch cell and is normalized over the
 * outgoing pitch grid. Reuse the rules across the table and the field
 * brackets across all columns on each flux surface.
 */
void MollerDeltaAngleKernel::TabulateDeltaMatrixOnXiStarGrid() {
    if (!hasActivePairs)
        return;

    const len_t Nr = gridK->GetNr();
    const KnockOnUtilities::Quadrature::GaussLegendreRule poloidalRule(poloidalOrder);
    const KnockOnUtilities::Quadrature::GaussLegendreRule pitchRule(pitchOrder);
    for (len_t ir = 0; ir < Nr; ++ir) {
        auto *fsa = gridP->GetRadialGrid()->GetFluxSurfaceAverager();
        KnockOnUtilities::Quadrature::OrbitGeometry geometry(fsa, ir);
        const len_t NxiK = gridK->GetNp2(ir);
        const len_t NxiP = gridP->GetNp2(ir);

        for (len_t m = 0; m < nXiStarsTabulate; ++m) {
            real_t *slice = deltaTable[ir][m];
            for (len_t l = 0; l < NxiP; ++l) {
                real_t *column = slice + l * NxiK;
                KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
                    ir, xiStarsTab[m], l, gridK, gridP, column, geometry, poloidalRule, pitchRule
                );
            }
        }
    }
}

// Precompute interpolation tables for linear interpolation in delta vs xi_star.
void MollerDeltaAngleKernel::BuildXiStarInterp() {
    const real_t inv_dXiStar = (dXiStar > 0) ? 1.0 / dXiStar : 0;
    const len_t n = nXiStarsTabulate;  // validated >= 2

    for (len_t ir = 0; ir < gridK->GetNr(); ++ir) {
        const auto *mgK = gridK->GetMomentumGrid(ir);
        const auto *mgP = gridP->GetMomentumGrid(ir);

        const len_t NpK = mgK->GetNp1();
        const len_t NpP = mgP->GetNp1();

        for (len_t i = 0; i < NpK; ++i) {
            for (len_t k = 0; k < NpP; ++k) {
                XiStarInterp &T = xiInterp[ir][i * NpP + k];
                T.m0 = 0;
                T.w1 = 0;
                real_t xs;
                if (!hasActivePairs || dXiStar == 0 || !EvaluateXiStarOnGrids(i, k, xs))
                    continue;

                if (xs <= xiStarMin) {
                    // clamp to first slice: D0=slice0, D1=slice1, w1=0 -> uses slice0 exactly
                    T.m0 = 0;
                    T.w1 = 0;
                    continue;
                }
                if (xs >= xiStarMax) {
                    // clamp to last slice: pick m0=n-2, w1=1 -> uses slice(n-1) exactly
                    T.m0 = n - 2;
                    T.w1 = 1;
                    continue;
                }

                const real_t s = (xs - xiStarMin) * inv_dXiStar;  // in (0, n-1)
                len_t m0 = (len_t)s;                              // floor
                if (m0 >= n - 1) {                                // safety
                    m0 = n - 2;
                }
                T.m0 = m0;
                T.w1 = s - (real_t)m0;
            }
        }
    }
}

void MollerDeltaAngleKernel::GridRebuilt() {
    Deallocate();
    ValidateGridAssumptions();

    BuildXiStarTabulationGrid();
    AllocateDeltaTables();
    TabulateDeltaMatrixOnXiStarGrid();
    BuildXiStarInterp();
}

/** Apply the angular kernel to all primary momenta for outgoing cell i. */
void MollerDeltaAngleKernel::AccumulatePitchAllPrimaries(
    len_t ir, len_t i, const real_t *W, const real_t *sigma,
    real_t *scratch, real_t *Cj
) const {
    if (!hasActivePairs) return;
    const len_t npP=gridP->GetNp1(ir);
    AccumulateMollerPitch(npP,gridP->GetNp2(ir),gridK->GetNp2(ir),
        nXiStarsTabulate,xiInterp[ir]+i*npP,sigma,W,deltaTable[ir],scratch,Cj);
}

/** Apply the same contraction to a single incident momentum cell. */
void MollerDeltaAngleKernel::AccumulatePitch(
    len_t ir, len_t i, len_t k, const real_t *W, real_t sigma, real_t *Cj
) const {
    if (sigma==0 || !hasActivePairs) return;
    const len_t nxiP=gridP->GetNp2(ir);
    const auto& t=xiInterp[ir][i*gridP->GetNp1(ir)+k];
    const XiStarInterp local{0,t.w1};
    std::vector<real_t> scratch(2*nxiP);
    AccumulateMollerPitch(1,nxiP,gridK->GetNp2(ir),2,&local,&sigma,W,
        deltaTable[ir]+t.m0,scratch.data(),Cj);
}
