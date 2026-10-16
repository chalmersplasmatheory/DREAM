/**
 * Implementation of AvalancheSourceDirect — direct numerical
 * integration of the Møller cross-section avalanche kernel.
 *
 * The kernel is pre-computed as a PETSc sparse matrix for each
 * radial grid point, and applied at runtime via sparse matrix-vector
 * multiplication.
 *
 * Physics reference:
 *   Rosenbluth & Putvinski, NF 37 (1997) 10
 *   Boozer, PPCF 57 (2015) 045001, SupplementMaterial
 *   Hesslow et al., NF 59 (2019) 084004
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceDirect.hpp"
#include "DREAM/Constants.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/MomentumGrid.hpp"
#include "FVM/Grid/RadialGrid.hpp"

#include <petscmat.h>
#include <cmath>

using namespace DREAM;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
AvalancheSourceDirect::AvalancheSourceDirect(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t scaleFactor
) : FVM::EquationTerm(kineticGrid),
    unknowns(u),
    pCutoff(pCutoff),
    scaleFactor(scaleFactor),
    Np1(0), Np2(0), NCellsPerRadius(0),
    nnz_per_row(nullptr),
    matricesBuilt(false),
    p_edges(nullptr),
    xi0_edges(nullptr)
{
    SetName("AvalancheSourceDirect");

    // Resolve unknown quantity IDs
    id_f_re   = u->GetUnknownID(OptionConstants::UQTY_F_RE);
    id_ntot   = u->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = u->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    // Add Jacobian dependencies
    AddUnknownForJacobian(u, id_ntot);
    AddUnknownForJacobian(u, id_Efield);

    // Pre-factor: e^4 / ((4πε₀)² m_e² c³) = r₀² c
    real_t e = Constants::ec;
    real_t epsmc = 4 * M_PI * Constants::eps0 * Constants::me * Constants::c;
    this->preFactor = (e * e * e * e) / (epsmc * epsmc * Constants::c);

    // Np1/Np2 may be zero if grid not yet initialized; GridRebuilt
    // will be called again by the framework when the grid is ready.
    this->GridRebuilt();
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
AvalancheSourceDirect::~AvalancheSourceDirect() {
    DeallocateMatrices();
    delete [] nnz_per_row;
}

//------------------------------------------------------------------------------
// Deallocate PETSc matrices
//------------------------------------------------------------------------------
void AvalancheSourceDirect::DeallocateMatrices() {
    // Guard against calling MatDestroy after PetscFinalize()
    // (destructor may run after dream_finalize() in main.cpp)
    PetscBool petscInit;
    PetscInitialized(&petscInit);
    if (!petscInit) {
        precomputedK.clear();
        allocated_nr = 0;
        matricesBuilt = false;
        return;
    }

    for (auto &mat : precomputedK) {
        if (mat != nullptr) {
            MatDestroy(&mat);
        }
    }
    precomputedK.clear();
    allocated_nr = 0;
    matricesBuilt = false;
}

//------------------------------------------------------------------------------
// Allocate PETSc matrices (one per radial position)
//------------------------------------------------------------------------------
void AvalancheSourceDirect::AllocateMatrices() {
    DeallocateMatrices();

    len_t nr = grid->GetNr();
    precomputedK.resize(nr, nullptr);

    // Momentum grid dimensions (assumed same at all radii)
    Np1 = grid->GetMomentumGrid(0)->GetNp1();
    Np2 = grid->GetMomentumGrid(0)->GetNp2();
    NCellsPerRadius = Np1 * Np2;

    // Cache grid edges
    p_edges   = grid->GetMomentumGrid(0)->GetP1_f();
    xi0_edges = grid->GetMomentumGrid(0)->GetP2_f();

    // Allocate/clear non-zero per row array
    delete [] nnz_per_row;
    nnz_per_row = new PetscInt[NCellsPerRadius];

    allocated_nr = nr;
}

//------------------------------------------------------------------------------
// Møller cross-section smooth factor A(p1, p)
//   Returns ntot * v_e * dσ/dp * (2π/c) * preFactor
//------------------------------------------------------------------------------
real_t AvalancheSourceDirect::ComputeA(real_t p1, real_t p, real_t ntot) const {
    real_t g1 = sqrt(1.0 + p1 * p1);
    real_t g  = sqrt(1.0 + p  * p);

    // Energy threshold: must satisfy 1 < γ < 2γ₁-1
    if (g - 1.0 < 1e-12 || g1 - 1.0 < 1e-12)
        return 0.0;
    if (g >= 2.0 * g1 - 1.0)
        return 0.0;

    real_t nu = (g - 1.0) / (g1 - 1.0);
    if (nu <= 0.0 || nu >= 1.0)
        return 0.0;

    real_t x = 1.0 / (nu * (1.0 - nu));

    // dσ/dν from Møller cross-section (relativistic)
    real_t dSigma_dNu = (g1 * g1) / ((g1 - 1.0) * (g1 - 1.0) * (g1 + 1.0))
                      * (x * x - 3.0 * x
                          + ((g1 - 1.0) * (g1 - 1.0) / (g1 * g1)) * (1.0 + x));

    // dν/dp
    real_t dNu_dp = p / (g * (g1 - 1.0));

    // dσ/dp
    real_t dSigma_dp = dSigma_dNu * dNu_dp;

    // v_e * dσ/dp * (2π/c) * preFactor
    real_t ve = p1 / g1;
    return ntot * ve * dSigma_dp * (2.0 * M_PI / Constants::c) * preFactor;
}

//------------------------------------------------------------------------------
// Compute angular support parameters ξ₁c and ξ₂
//   Using the Boozer (2015) SupplementMaterial formula.
//------------------------------------------------------------------------------
void AvalancheSourceDirect::ComputeAngularSupport(
    real_t p1, real_t xi01, real_t p,
    real_t &xi1c, real_t &xi2
) {
    real_t g1 = sqrt(1.0 + p1 * p1);
    real_t g  = sqrt(1.0 + p  * p);

    real_t gamma_rel = (g1 + 1.0) * (g - 1.0) / ((g1 - 1.0) * (g + 1.0));
    if (gamma_rel < 0.0) {
        xi1c = 0.0;
        xi2  = -1.0;   // marker for invalid
        return;
    }
    real_t xi_star = sqrt(gamma_rel);

    // Center and width of the angular distribution
    xi1c = xi01 * xi_star;
    xi2  = sqrt(fmax(0.0, 1.0 - xi01 * xi01))
         * sqrt(fmax(0.0, 1.0 - xi_star * xi_star));

    if (xi2 < 1e-12)
        xi2 = 1e-12;
}

//------------------------------------------------------------------------------
// Static callback for EvaluateCellAveragedBounceIntegralOverP2.
//
// Returns the value of the angular kernel Π at a given poloidal angle θ.
//------------------------------------------------------------------------------
real_t AvalancheSourceDirect::F_AngularKernel(
    real_t xiOverXi0, real_t BOverBmin,
    real_t ROverR0, real_t NablaR2, void *par
) {
    (void)ROverR0;
    (void)NablaR2;
    (void)BOverBmin;

    auto *params = static_cast<AngularKernelParams*>(par);

    // Recover local pitch ξ at this poloidal angle
    real_t xi0_val = params->xi0;
    real_t xi  = xiOverXi0 * xi0_val;

    // Source pitch: use bounce-invariant pitch directly
    real_t xi1 = params->xi1_c;

    // Compute the angular support
    real_t xi1c, xi2;
    ComputeAngularSupport(params->p1, xi1, params->p, xi1c, xi2);

    if (xi2 < 0.0)
        return 0.0;

    // Π function: 1/(π √(ξ₂² - (ξ - ξ₁c)²))
    real_t dxi = xi - xi1c;
    real_t arg = xi2 * xi2 - dxi * dxi;
    if (arg <= 0.0)
        return 0.0;

    return 1.0 / (M_PI * sqrt(arg));
}

//------------------------------------------------------------------------------
// Compute bounce-averaged angular kernel for a target ξ-cell [xi_l, xi_u]
// and source ξ-cell centered at xi1_c.
//
// Returns the cell-averaged bounce integral divided by Vp, giving
// the flux-surface-averaged angular kernel per unit ξ₀.
//------------------------------------------------------------------------------
real_t AvalancheSourceDirect::ComputeAngularKernel(
    len_t ir, real_t p,
    real_t xi_l, real_t xi_u,
    real_t p1, real_t xi1_c, real_t xi1_l, real_t xi1_u,
    int_t RESign
) {
    (void)RESign;

    if (fabs(xi_u - xi_l) < 100.0 * std::numeric_limits<real_t>::epsilon())
        return 0.0;

    real_t xi0_c = 0.5 * (xi_l + xi_u);

    DREAM::FVM::FluxSurfaceAverager *fsa = this->GetFSA(ir);
    if (fsa == nullptr)
        return 0.0;

    // Prepare parameters for the callback
    AngularKernelParams angParams;
    angParams.p     = p;
    angParams.p1    = p1;
    angParams.xi0   = xi0_c;
    angParams.xi1_c = xi1_c;
    angParams.dxi1  = xi1_u - xi1_l;
    angParams.xi1_l = xi1_l;
    angParams.xi1_u = xi1_u;
    angParams.ir    = ir;
    angParams.fsa   = this->GetFSA(ir);
    angParams.RESign = RESign;

    // Evaluate the bounce integral averaged over the target ξ-cell.
    // Returns ⟨∮ F √g dθ⟩_cell / dxi_cell
    real_t cellAvgBounceIntegral = fsa->EvaluateCellAveragedBounceIntegralOverP2(
        ir, xi_l, xi_u, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION,
        &F_AngularKernel, &angParams
    );

    // Find the Vp at this xi cell for normalization
    const real_t *xi0_f2_local = grid->GetMomentumGrid(ir)->GetP2_f();
    len_t j_xi = 0;
    for (len_t j = 0; j < Np2; j++) {
        if (xi0_c >= xi0_f2_local[j] && xi0_c <= xi0_f2_local[j+1]) {
            j_xi = j;
            break;
        }
    }
    // Vp at cell center (p=0 index, independent of p in p-xi grid)
    real_t Vp_val_local = grid->GetVp(ir, 0, j_xi);

    // The bounce integral / Vp gives the bounce average
    if (Vp_val_local <= 0.0)
        return 0.0;

    return cellAvgBounceIntegral / Vp_val_local;
}

//------------------------------------------------------------------------------
// Build the kernel matrix for a single radial index.
//
// K[ir] has dimensions NCellsPerRadius × NCellsPerRadius.
// Row = target cell (jt * Np1 + it), Col = source cell (js * Np1 + is).
//------------------------------------------------------------------------------
void AvalancheSourceDirect::BuildKernelMatrixForRadialIndex(len_t ir) {
    auto *mg = grid->GetMomentumGrid(ir);
    const real_t *p_edges_ir   = mg->GetP1_f();
    const real_t *xi0_edges_ir = mg->GetP2_f();
    const real_t *xi0_cells    = mg->GetP2();

    len_t np1 = mg->GetNp1();
    len_t np2 = mg->GetNp2();
    len_t nCells = np1 * np2;

    // Compute trapped-passing boundary: only passing particles (xi > xi_trapped)
    // contribute to avalanche generation
    const real_t *Bmin_arr = grid->GetRadialGrid()->GetBmin();
    const real_t *Bmax_arr = grid->GetRadialGrid()->GetBmax();
    real_t xi_trapped = sqrt(fmax(0.0, 1.0 - Bmin_arr[ir] / fmax(Bmax_arr[ir], 1e-30)));

    // Phase 1: count non-zeros per row
    PetscInt *d_nnz = new PetscInt[nCells]();

    for (len_t it = 0; it < np1; it++) {
        real_t p_mid = 0.5 * (p_edges_ir[it] + p_edges_ir[it+1]);

        for (len_t jt = 0; jt < np2; jt++) {
            // Skip target cells in trapped region
            if (xi0_edges_ir[jt+1] <= xi_trapped)
                continue;

            len_t row = jt * np1 + it;

            for (len_t is = 0; is < np1; is++) {
                real_t p1_mid = 0.5 * (p_edges_ir[is] + p_edges_ir[is+1]);
                if (p1_mid < pCutoff)
                    continue;

                for (len_t js = 0; js < np2; js++) {
                    real_t xi01 = xi0_cells[js];
                    // Skip source cells in trapped region
                    if (xi01 <= xi_trapped)
                        continue;

                    real_t xi1c, xi2;
                    ComputeAngularSupport(p1_mid, xi01, p_mid, xi1c, xi2);
                    if (xi2 < 0.0)
                        continue;

                    // Skip if source Vp is zero (boundary cells)
                    if (grid->GetVp(ir, is, js) <= 0.0)
                        continue;

                    // Check overlap with target xi-cell
                    real_t xi_l = xi0_edges_ir[jt];
                    real_t xi_u = xi0_edges_ir[jt+1];
                    real_t a = fmax(xi_l, xi1c - xi2);
                    real_t b = fmin(xi_u, xi1c + xi2);
                    if (b <= a)
                        continue;

                    d_nnz[row]++;
                }
            }
        }
    }

    // Create PETSc matrix with pre-allocation
    Mat K;
    MatCreateSeqAIJ(PETSC_COMM_SELF, (PetscInt)nCells, (PetscInt)nCells,
                    0, d_nnz, &K);
    MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // Phase 2: fill matrix elements
    for (len_t it = 0; it < np1; it++) {
        real_t p_cen  = 0.5 * (p_edges_ir[it] + p_edges_ir[it+1]);
        real_t dp_t   = p_edges_ir[it+1] - p_edges_ir[it];

        for (len_t jt = 0; jt < np2; jt++) {
            // Skip target cells in trapped region
            if (xi0_edges_ir[jt+1] <= xi_trapped)
                continue;

            len_t row = jt * np1 + it;
            real_t xi_l = xi0_edges_ir[jt];
            real_t xi_u = xi0_edges_ir[jt+1];
            real_t dxi_t = xi_u - xi_l;

            real_t Vp_target = grid->GetVp(ir, it, jt);

            for (len_t is = 0; is < np1; is++) {
                real_t p1_cen = 0.5 * (p_edges_ir[is] + p_edges_ir[is+1]);
                real_t dp1_s  = p_edges_ir[is+1] - p_edges_ir[is];

                // Energy threshold
                real_t g1 = sqrt(1.0 + p1_cen * p1_cen);
                real_t g  = sqrt(1.0 + p_cen  * p_cen);
                if (g >= 2.0 * g1 - 1.0 || g >= g1 - 1e-12)
                    continue;
                if (p1_cen < pCutoff)
                    continue;

                for (len_t js = 0; js < np2; js++) {
                    real_t xi01 = xi0_cells[js];
                    // Skip source cells in trapped region
                    if (xi01 <= xi_trapped)
                        continue;

                    len_t col = js * np1 + is;
                    real_t xi1_l = xi0_edges_ir[js];
                    real_t xi1_u = xi0_edges_ir[js+1];
                    real_t xi1_c = 0.5 * (xi1_l + xi1_u);
                    real_t dxi_s = xi1_u - xi1_l;

                    // Compute angular overlap
                    real_t xi1c_ang, xi2_ang;
                    ComputeAngularSupport(p1_cen, xi1_c, p_cen,
                                          xi1c_ang, xi2_ang);
                    if (xi2_ang < 0.0)
                        continue;

                    real_t a = fmax(xi_l, xi1c_ang - xi2_ang);
                    real_t b = fmin(xi_u, xi1c_ang + xi2_ang);
                    if (b <= a)
                        continue;

                    // Analytical integral of Π over target xi-cell
                    // Clamp asin arguments to [-1, 1] to protect against
                    // floating-point rounding that would produce NaN.
                    real_t arg_a = (a - xi1c_ang) / xi2_ang;
                    real_t arg_b = (b - xi1c_ang) / xi2_ang;
                    arg_a = fmax(-1.0, fmin(1.0, arg_a));
                    arg_b = fmax(-1.0, fmin(1.0, arg_b));
                    real_t u_a = asin(arg_a);
                    real_t u_b = asin(arg_b);
                    real_t xi_int = (u_b - u_a) / M_PI;
                    if (xi_int <= 0.0)
                        continue;

                    // Møller factor (ntot=1, scaled at runtime)
                    real_t A_val = ComputeA(p1_cen, p_cen, 1.0);
                    if (A_val == 0.0)
                        continue;

                    // Source cell volume in bounce-averaged phase space
                    real_t Vp_source = grid->GetVp(ir, is, js);
                    if (Vp_source <= 0.0)
                        continue;

                    // Matrix element with volume normalization
                    real_t vol_source = Vp_source * p1_cen * p1_cen * dp1_s * dxi_s;
                    real_t vol_target = Vp_target * p_cen * p_cen * dp_t * dxi_t;
                    // Skip if target volume is too small relative to source
                    // to avoid producing Inf/unstable matrix entries.
                    if (vol_target <= 1e-30 * vol_source)
                        continue;

                    real_t val = A_val * xi_int * vol_source / vol_target;

                    // Guard against NaN/Inf from any edge case
                    if (!std::isfinite(val))
                        continue;

                    PetscInt row_p = (PetscInt)row;
                    PetscInt col_p = (PetscInt)col;
                    PetscScalar val_p = (PetscScalar)val;
                    MatSetValues(K, 1, &row_p, 1, &col_p, &val_p, INSERT_VALUES);
                }
            }
        }
    }

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

    precomputedK[ir] = K;
    delete [] d_nnz;
}

//------------------------------------------------------------------------------
// Build all kernel matrices
//------------------------------------------------------------------------------
void AvalancheSourceDirect::BuildKernelMatrices() {
    len_t nr = grid->GetNr();
    for (len_t ir = 0; ir < nr; ir++) {
        BuildKernelMatrixForRadialIndex(ir);
    }
    matricesBuilt = true;
}

//------------------------------------------------------------------------------
// GridRebuilt — called when the computational grid changes
//------------------------------------------------------------------------------
bool AvalancheSourceDirect::GridRebuilt() {
    FVM::EquationTerm::GridRebuilt();
    AllocateMatrices();
    BuildKernelMatrices();
    return true;
}

//------------------------------------------------------------------------------
// Rebuild — called at the beginning of each time step
//------------------------------------------------------------------------------
void AvalancheSourceDirect::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {
    // Kernel matrices are pre-computed in GridRebuilt().
    // Runtime scaling by n_tot is done in SetVectorElements.
}

//------------------------------------------------------------------------------
// SetVectorElements — evaluate the source term explicitly:
//   vec += S = n_tot * scaleFactor * K_matrix * f_re
//------------------------------------------------------------------------------
void AvalancheSourceDirect::SetVectorElements(
    real_t *vec, const real_t *x
) {
    if (!matricesBuilt)
        return;

    const real_t *f_re = this->unknowns->GetUnknownData(id_f_re);
    const real_t *ntot = this->unknowns->GetUnknownData(id_ntot);

    Vec f_vec, S_vec;
    VecCreateSeq(PETSC_COMM_SELF, (PetscInt)NCellsPerRadius, &f_vec);
    VecCreateSeq(PETSC_COMM_SELF, (PetscInt)NCellsPerRadius, &S_vec);

    len_t offset = 0;
    for (len_t ir = 0; ir < static_cast<len_t>(precomputedK.size()); ir++) {
        if (precomputedK[ir] == nullptr) {
            offset += NCellsPerRadius;
            continue;
        }

        // Copy current iterate into PETSc vector
        const real_t *src = (x != nullptr) ? (x + offset) : (f_re + offset);
        PetscScalar *f_array;
        VecGetArray(f_vec, &f_array);
        for (len_t i = 0; i < NCellsPerRadius; i++)
            f_array[i] = (PetscScalar)src[i];
        VecRestoreArray(f_vec, &f_array);

        // S = K * f
        MatMult(precomputedK[ir], f_vec, S_vec);

        // Add to global vector: vec += n_tot * scaleFactor * S
        PetscScalar *s_array;
        VecGetArray(S_vec, &s_array);
        real_t ntot_ir = ntot[ir];
        for (len_t i = 0; i < NCellsPerRadius; i++)
            vec[offset + i] += (real_t)(scaleFactor * ntot_ir * s_array[i]);
        VecRestoreArray(S_vec, &s_array);

        offset += NCellsPerRadius;
    }

    VecDestroy(&f_vec);
    VecDestroy(&S_vec);
}

//------------------------------------------------------------------------------
// SetMatrixElements — semi-implicit matrix:
//   mat += n_tot * scaleFactor * K_matrix
//------------------------------------------------------------------------------
void AvalancheSourceDirect::SetMatrixElements(
    FVM::Matrix *mat, real_t * /*rhs*/
) {
    if (!matricesBuilt)
        return;

    const real_t *ntot = this->unknowns->GetUnknownData(id_ntot);

    len_t offset = 0;
    for (len_t ir = 0; ir < static_cast<len_t>(precomputedK.size()); ir++) {
        if (precomputedK[ir] == nullptr) {
            offset += NCellsPerRadius;
            continue;
        }

        real_t ntot_ir = ntot[ir];

        for (len_t row_local = 0; row_local < NCellsPerRadius; row_local++) {
            PetscInt ncols;
            const PetscInt *cols;
            const PetscScalar *vals;
            MatGetRow(precomputedK[ir], (PetscInt)row_local,
                      &ncols, &cols, &vals);

            if (ncols > 0) {
                std::vector<PetscInt> local_cols(cols, cols + ncols);
                std::vector<PetscScalar> scaled_vals(ncols);
                for (PetscInt c = 0; c < ncols; c++)
                    scaled_vals[c] = (PetscScalar)(scaleFactor * ntot_ir * vals[c]);

                mat->SetRow((PetscInt)row_local, ncols, local_cols.data(),
                            scaled_vals.data(), ADD_VALUES);
            }

            MatRestoreRow(precomputedK[ir], (PetscInt)row_local,
                          &ncols, &cols, &vals);
        }

        offset += NCellsPerRadius;
    }
}

//------------------------------------------------------------------------------
// SetJacobianBlock — provide analytic Jacobian for the Newton solver.
//------------------------------------------------------------------------------
bool AvalancheSourceDirect::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId,
    FVM::Matrix *jac, const real_t *x
) {
    if (!matricesBuilt)
        return false;

    // f_re-f_re block: the kernel matrix itself
    if (uqtyId == id_f_re && derivId == id_f_re) {
        SetMatrixElements(jac, nullptr);
        return true;
    }

    // n_tot derivatives: diagonal contributions in radius
    if (derivId == id_ntot) {
        if (x == nullptr) return false;

        const real_t *f_re = this->unknowns->GetUnknownData(id_f_re);

        Vec f_vec, S_vec;
        VecCreateSeq(PETSC_COMM_SELF, (PetscInt)NCellsPerRadius, &f_vec);
        VecCreateSeq(PETSC_COMM_SELF, (PetscInt)NCellsPerRadius, &S_vec);

        len_t offset = 0;
        for (len_t ir = 0; ir < static_cast<len_t>(precomputedK.size()); ir++) {
            if (precomputedK[ir] == nullptr) {
                offset += NCellsPerRadius;
                continue;
            }

            PetscScalar *f_array;
            VecGetArray(f_vec, &f_array);
            for (len_t i = 0; i < NCellsPerRadius; i++)
                f_array[i] = (PetscScalar)f_re[offset + i];
            VecRestoreArray(f_vec, &f_array);

            MatMult(precomputedK[ir], f_vec, S_vec);

            PetscScalar *s_array;
            VecGetArray(S_vec, &s_array);
            for (len_t i = 0; i < NCellsPerRadius; i++) {
                PetscInt row = (PetscInt)i;
                PetscInt col = (PetscInt)ir;
                PetscScalar val = (PetscScalar)(scaleFactor * s_array[i]);
                jac->SetElement(row, col, val);
            }
            VecRestoreArray(S_vec, &s_array);

            offset += NCellsPerRadius;
        }

        VecDestroy(&f_vec);
        VecDestroy(&S_vec);

        return true;
    }

    // E-field derivatives
    if (derivId == id_Efield)
        return false;

    return false;
}
