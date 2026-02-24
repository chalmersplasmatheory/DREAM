#include "DREAM/Equations/Kinetic/MollerBoltzmannOperator.hpp"

#include <cmath>
#include <limits>

#include "DREAM/DREAMException.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

using namespace DREAM;

MollerBoltzmannOperator::MollerBoltzmannOperator(
    FVM::Grid *gridKnockon, const FVM::Grid *grid_primary, FVM::UnknownQuantityHandler *unknowns,
    len_t id_f_primary, const MollerEnergyKernel *energyKernel,
    const MollerDeltaAngleKernel *angleKernel, real_t scaleFactor
)
    : FVM::EquationTerm(gridKnockon),
      gridPrimary(grid_primary),
      unknowns(unknowns),
      id_f_primary(id_f_primary),
      scaleFactor(scaleFactor),
      energyKernel(energyKernel),
      angleKernel(angleKernel) {
    SetName("MollerBoltzmannOperator");

    ValidateInputParameters();
    ValidateGridAssumptions();

    const len_t id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    AddUnknownForJacobian(unknowns, id_ntot);

    t_source_rebuilt = -std::numeric_limits<real_t>::infinity();

    AllocateSourceVector();
    AllocateScratchBuffers();
}

MollerBoltzmannOperator::~MollerBoltzmannOperator() { Deallocate(); }

void MollerBoltzmannOperator::ValidateInputParameters() const {
    if (unknowns == nullptr) {
        throw DREAMException("MollerBoltzmannOperator: 'unknowns' must not be null.");
    }
}

void MollerBoltzmannOperator::ValidateGridAssumptions() const {
    if (grid == nullptr || gridPrimary == nullptr) {
        throw DREAMException("MollerBoltzmannOperator: grid pointers must not be null.");
    }
    if (grid->GetNr() != gridPrimary->GetNr()) {
        throw DREAMException("MollerBoltzmannOperator: grid and grid_primary must have same Nr.");
    }
    for (len_t ir = 1; ir < grid->GetNr(); ++ir) {
        if (grid->GetMomentumGrid(ir) != grid->GetMomentumGrid(0)) {
            throw DREAMException(
                "MollerBoltzmannOperator: requires shared knock-on MomentumGrid across radii."
            );
        }
        if (gridPrimary->GetMomentumGrid(ir) != gridPrimary->GetMomentumGrid(0)) {
            throw DREAMException(
                "MollerBoltzmannOperator: requires shared primary MomentumGrid across radii."
            );
        }
    }
}

void MollerBoltzmannOperator::AllocateSourceVector() {
    sourceVector = new real_t[grid->GetNCells()];
}

void MollerBoltzmannOperator::AllocateScratchBuffers() {
    const auto *mgK = grid->GetMomentumGrid(0);
    const auto *mgP = gridPrimary->GetMomentumGrid(0);

    const len_t Np1P = mgP->GetNp1();
    const len_t NxiP = mgP->GetNp2();
    const len_t NxiK = mgK->GetNp2();

    primaryWeights = new real_t[Np1P * NxiP];
    Cj = new real_t[NxiK];
}

void MollerBoltzmannOperator::Deallocate() {
    if (sourceVector != nullptr) {
        delete[] sourceVector;
        sourceVector = nullptr;
    }
    if (primaryWeights != nullptr) {
        delete[] primaryWeights;
        primaryWeights = nullptr;
    }
    if (Cj != nullptr) {
        delete[] Cj;
        Cj = nullptr;
    }
}


// Optimization: since the full equation term is given by
//   S_ij = \int dp1 \int dxi01 S_i * delta_jk * Vp * f
//        ~ sum_k sum_l dp_k dxi_l S_ik delta_jklm Vp_kl f_kl
//        = sum_k sum_l W_kl S_ik delta_jklm
// we here collect the purely "p1"-dependent weight 
//   W_kl = dp_k * dxi_l * Vp_kl * f_kl 
// This saves us from recomputing this product for each (i, j).
void MollerBoltzmannOperator::BuildPrimaryWeights(
    len_t ir, const real_t *f_primary_ir, real_t *W_k_l
) const {
    const auto *mgP = gridPrimary->GetMomentumGrid(ir);

    const len_t Np1P = mgP->GetNp1();
    const len_t NxiP = mgP->GetNp2();

    const real_t *dp1 = mgP->GetDp1();
    const real_t *dxi1 = mgP->GetDp2();
    const real_t *VpP = gridPrimary->GetVp(ir);

    for (len_t l = 0; l < NxiP; ++l) {
        const real_t dxi = dxi1[l];
        const len_t base_lk = l * Np1P;
        for (len_t k = 0; k < Np1P; ++k) {
            const len_t idx_lk = base_lk + k;
            W_k_l[k * NxiP + l] = dp1[k] * dxi * VpP[idx_lk] * f_primary_ir[idx_lk];
        }
    }
}

/**
 * Assemble the explicit-time Boltzmann source on the knock-on grid from f_primary.
 *
 * Key idea:
 *  - First precompute W_{k l} = dp1_k dxi_l Vp_{k l} f_{k l} for each radius,
 *    then for each outgoing momentum cell i, accumulate C_j by looping over k and
 *    applying the interpolated delta kernel over pitch.
 */
void MollerBoltzmannOperator::SetSourceVector(const real_t *f_primary) {
    len_t offK = 0;
    len_t offP = 0;

    const len_t Nr = grid->GetNr();

    for (len_t ir = 0; ir < Nr; ++ir) {
        const auto *mgK = grid->GetMomentumGrid(ir);
        const auto *mgP = gridPrimary->GetMomentumGrid(ir);

        const len_t Np1K = mgK->GetNp1();
        const len_t NxiK = mgK->GetNp2();
        const len_t Np1P = mgP->GetNp1();
        const len_t NxiP = mgP->GetNp2();

        const real_t *VpK = grid->GetVp(ir);

        // W_{k,l} = dp1_k * dxi_l * Vp_{k,l} * f_{k,l}
        BuildPrimaryWeights(ir, f_primary + offP, primaryWeights);
        for (len_t i = 0; i < Np1K; ++i) {
            for (len_t j = 0; j < NxiK; ++j) {
                Cj[j] = 0.0;
            }
            for (len_t k = 0; k < Np1P; ++k) {
                const real_t sigma_ik = energyKernel->DifferentialCS(i, k);
                const real_t *Wkl = primaryWeights + k * NxiP;
                // for each j (outgoing), let the angleKernel carry out
                // the pitch integration over l (incident)
                angleKernel->AccumulatePitch(ir, i, k, Wkl, sigma_ik, Cj);
            }
            for (len_t j = 0; j < NxiK; ++j) {
                const len_t idxK = offK + j * Np1K + i;
                const real_t Vp = VpK[j * Np1K + i];
                sourceVector[idxK] = (Vp != 0) ? (Cj[j] * (scaleFactor / Vp)) : 0.0;
            }
        }
        offK += mgK->GetNCells();
        offP += mgP->GetNCells();
    }
}

void MollerBoltzmannOperator::Rebuild(real_t t, real_t /*dt*/, FVM::UnknownQuantityHandler *uqh) {
    if (uqh != nullptr) {
        this->unknowns = uqh;
    }
    constexpr real_t TIME_EPS_FACTOR = 100;
    const bool timeHasUpdated =
        std::abs(t - t_source_rebuilt) > TIME_EPS_FACTOR * std::numeric_limits<real_t>::epsilon();

    // only update the f_primary-dependent source function once per time step
    if (timeHasUpdated) {
        const real_t *f_primary = unknowns->GetUnknownDataPrevious(id_f_primary);
        SetSourceVector(f_primary);
        t_source_rebuilt = t;
    }
}

void MollerBoltzmannOperator::SetVectorElements(real_t *vec, const real_t *x) {
    len_t offset = 0;
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        const auto *mg = grid->GetMomentumGrid(ir);
        const len_t Np = mg->GetNp1();
        const len_t Nxi = mg->GetNp2();

        for (len_t j = 0; j < Nxi; ++j) {
            const len_t base = offset + Np * j;
            for (len_t i = 0; i < Np; ++i) {
                const len_t ind = base + i;
                vec[ind] += x[ir] * sourceVector[ind];
            }
        }
        offset += mg->GetNCells();
    }
}

bool MollerBoltzmannOperator::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t * /*x*/
) {
    if (derivId != uqtyId) {
        return false;
    }
    len_t offset = 0;
    for (len_t ir = 0; ir < grid->GetNr(); ++ir) {
        const auto *mg = grid->GetMomentumGrid(ir);
        const len_t Np = mg->GetNp1();
        const len_t Nxi = mg->GetNp2();

        for (len_t j = 0; j < Nxi; ++j) {
            const len_t base = offset + Np * j;
            for (len_t i = 0; i < Np; ++i) {
                const len_t ind = base + i;
                jac->SetElement(ind, ir, sourceVector[ind]);
            }
        }
        offset += mg->GetNCells();
    }

    return true;
}

bool MollerBoltzmannOperator::GridRebuilt() {
    Deallocate();
    ValidateGridAssumptions();

    // Reset explicit rebuild time.
    t_source_rebuilt = -std::numeric_limits<real_t>::infinity();

    // Reallocate operator-owned storage sized to new grid.
    AllocateSourceVector();
    AllocateScratchBuffers();

    return true;
}
