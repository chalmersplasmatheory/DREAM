/**
 * Implementation of a kinetic avalanche source term which uses
 * direct numerical integration of the Møller cross-section kernel,
 * without relying on the Rosenbluth-Putvinski delta-function
 * approximation in pitch angle.
 *
 * The source term takes the form
 *     S(p,xi) = n_tot(r) * ∫∫ f_re(p1,xi1) * K(p1,xi1 → p,xi) * p1² dp1 dxi1
 * where K = A(p1,p) * Π(xi; p1,xi1,p) contains the Møller cross-section
 * and the angular distribution function.
 *
 * The kernel is pre-computed as a sparse matrix (one per radial grid point)
 * and applied at runtime via sparse matrix-vector multiplication.
 */

#ifndef _DREAM_EQUATIONS_AVALANCHE_SOURCE_DIRECT_HPP
#define _DREAM_EQUATIONS_AVALANCHE_SOURCE_DIRECT_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include "FVM/Grid/FluxSurfaceAverager.hpp"
#include <petscmat.h>
#include <limits>
#include <vector>

namespace DREAM {

    class AvalancheSourceDirect : public FVM::EquationTerm {
    private:
        // Unknown quantity handler (stored locally since we don't inherit from FluidSourceTerm)
        FVM::UnknownQuantityHandler *unknowns;

        // Unknown quantity IDs
        len_t id_f_re;
        len_t id_ntot;
        len_t id_Efield;

        // Physical parameters
        real_t pCutoff;
        real_t scaleFactor;
        real_t preFactor;   // e^4 / ((4πε₀)² m_e² c³)

        // Grid dimensions (assumed same at all radii, as is standard in DREAM)
        len_t Np1, Np2;
        len_t NCellsPerRadius;   // Np1 * Np2

        // Pre-computed kernel matrices — one PETSc Mat per radial grid point
        std::vector<Mat> precomputedK;
        len_t allocated_nr;      // number of radial positions with allocated matrices
        // Exact non-zero count per row [NCellsPerRadius]
        PetscInt *nnz_per_row;
        bool matricesBuilt;

        // Cached grid geometry
        const real_t *p_edges;    // [Np1+1]
        const real_t *xi0_edges;  // [Np2+1]

        // Bounce-averager accessor shortcut
        DREAM::FVM::FluxSurfaceAverager *GetFSA(len_t /*ir*/) const {
            return grid->GetRadialGrid()->GetFluxSurfaceAverager();
        }

        // ---------- Internal methods ----------

        /**
         * Build the pre-computed kernel matrices for all radii.
         * Called once when the grid is (re)built.
         */
        void BuildKernelMatrices();

        /**
         * Build the kernel matrix for a single radial index.
         */
        void BuildKernelMatrixForRadialIndex(len_t ir);

        /**
         * Møller cross-section smooth factor A(p1, p).
         * Returns the part of the kernel that depends only on momentum magnitudes.
         */
        real_t ComputeA(real_t p1, real_t p, real_t ntot) const;

        /**
         * Compute the bounce-averaged angular kernel integrated over
         * the target xi-cell [xi_l, xi_u] and averaged over the source
         * xi-cell [xi1_l, xi1_u].
         *
         * Uses EvaluateCellAveragedBounceIntegralOverP2 internally.
         */
        real_t ComputeAngularKernel(
            len_t ir, real_t p,
            real_t xi_l, real_t xi_u,
            real_t p1, real_t xi1_c, real_t xi1_l, real_t xi1_u,
            int_t RESign
        );

        /**
         * Parameters passed into F_AngularKernel callback.
         */
        struct AngularKernelParams {
            real_t p, p1;           // momentum magnitudes
            real_t xi0;             // bounce-invariant pitch of target cell
            real_t xi1_c;           // center of source xi-cell (in xi₀)
            real_t dxi1;            // width of source xi-cell
            real_t xi1_l, xi1_u;    // bounds of source xi-cell
            len_t ir;               // radial index
            DREAM::FVM::FluxSurfaceAverager *fsa;
            int_t RESign;           // sign of electric field
        };

        /**
         * Static callback for EvaluateCellAveragedBounceIntegralOverP2.
         * Returns the angular kernel value at a given poloidal angle.
         */
        static real_t F_AngularKernel(
            real_t xiOverXi0, real_t BOverBmin,
            real_t ROverR0, real_t NablaR2, void *par
        );

        /**
         * Compute the pitch-angle support parameters ξ₁c and ξ₂
         * from Boozer kinematics.
         *
         * ξ₁c = center of the pitch distribution of created electrons
         * ξ₂  = half-width of the pitch distribution
         */
        static void ComputeAngularSupport(
            real_t p1, real_t xi01, real_t p,
            real_t &xi1c, real_t &xi2
        );

        // PETSc matrix management
        void AllocateMatrices();
        void DeallocateMatrices();

    protected:
        // Number of non-zeros per row (used for PETSc pre-allocation)
        virtual len_t GetNumberOfNonZerosPerRow() const override
        { return NCellsPerRadius; } // upper bound, refined in constructor

    public:
        AvalancheSourceDirect(
            FVM::Grid *kineticGrid,
            FVM::UnknownQuantityHandler *u,
            real_t pCutoff,
            real_t scaleFactor
        );
        ~AvalancheSourceDirect();

        // Grid rebuilt
        virtual bool GridRebuilt() override;

        // Rebuild term at new time step
        virtual void Rebuild(
            const real_t, const real_t, FVM::UnknownQuantityHandler*
        ) override;

        // Matrix / vector interface
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual bool SetJacobianBlock(
            const len_t uqtyId, const len_t derivId,
            FVM::Matrix *jac, const real_t *x
        ) override;
    };

}

#endif /* _DREAM_EQUATIONS_AVALANCHE_SOURCE_DIRECT_HPP */
