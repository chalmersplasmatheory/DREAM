#ifndef _DREAM_FVM_ADVECTION_INTERPOLATION_COEFFICIENT_TERM_HPP
#define _DREAM_FVM_ADVECTION_INTERPOLATION_COEFFICIENT_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

namespace DREAM::FVM {
    class AdvectionInterpolationCoefficient{
    public:

/**
 * Available interpolation methods. The expressions below are 
 * for constant grid spacing, but are implemented for general grids. 
 * Indices are relative to the flow; for positive advection coefficients
 *       y_{i-k} = y_[i-k]
 *       whereas for negatve advection coefficients
 *       y_{i-k} = y[i+k-1].
 *  CENTRED: Central difference scheme (2nd order, prone to oscillation) 
 *              y_{i-1/2} = (y_{i-1} + y_i)/2
 *  UPWIND: First order upwind (first order, unconditionally stable)
 *              y_{i-1/2} = y_{i-1}
 *  UPWIND_2ND_ORDER: Second order upwind (second order, but not bounded)
 *              y_{i-1/2} = (3*y_{i-1} - y_{i-2})/2
 *  DOWNWIND: First order downwind (first order, typically unphysical since non-transportive but can support strong gradients)
 *              y_{i-1/2} = y_i
 *  QUICK: Quadratic upwind (third order, relatively oscillation-resistant but not bounded)
 *              y_{i-1/2} = (3/4)y_{i-1} + (3/8)y_{i-2} - (1/8)y_i
 *  SMART: A flux-limited scheme that maximizes the use of the accurate QUICK scheme,
 *         but applies upwind in regions of strong gradients (where QUICK overshoots) 
 *         and matches to another scheme in a transition region in order to make the 
 *         interpolation coefficients continuous functions of flux and solution field.
 *         Benefit:  Bounded (essentially positivity-conserving), yet accurate since it 
 *                   maximizes the use of the third-order QUICK scheme. 
 *         Drawback: Discontinuously switches between interpolation schemes, which can
 *                   require additional iterations before convergence is reached (or in 
 *                   principle even not converge, where oscillation between two solutions can occur).
 *         Note that our implementation of the flux limited scheme for non-uniform grids
 *         differs from common expositions in the literature; while the method implemented in
 *         DREAM will not be strictly bounded, the QUICK scheme that is used in a majority
 *         of the domain will be correctly third-order accurate (unlike commonly cited schemes).
 *         In regions of constant grid spacing, the implemented scheme is exactly original SMART.
 *  MUSCL: A flux limited scheme, where the flux-limiter function takes the form
 *              psi(r) = max( 0, min(2*r, 0.5+0.5*r, 2) )
 *         Typically, the method is expected to be slightly less accurate than SMART but converges faster.
 *  OSPRE: A continuous flux limiter given by 1.5*r*(r+1)/(r*r+r+1). Converges much more rapidly
 *         than the piecewise linear limiters.
 *  TCDF:  A piecewise but continuously differentiable limiter, which corresponds to QUICK
 *         in a large part of the solution domain. Seems like the best of both worlds?
 */
        enum adv_interpolation {
            AD_INTERP_CENTRED  = 1,
            AD_INTERP_UPWIND   = 2,
            AD_INTERP_UPWIND_2ND_ORDER = 3,
            AD_INTERP_DOWNWIND = 4,
            AD_INTERP_QUICK = 5,
            AD_INTERP_SMART = 6,
            AD_INTERP_MUSCL = 7,
            AD_INTERP_OSPRE = 8,
            AD_INTERP_TCDF  = 9
        };


        enum adv_interp_mode {
            AD_INTERP_MODE_FULL,
            AD_INTERP_MODE_JACOBIAN
        };

        /** 
         * Default boundary conditions that can be used for
         * advection terms: Mirrored means that 
         *  f(x_min-x) = f(x_min + x) 
         * or 
         *  f(x_max - x) = f(x_max + x), 
         * whereas dirichlet is 
         *  f=0 
         * on boundary (conservative; zero flux)
         */
        enum adv_bc {
            AD_BC_MIRRORED,
            AD_BC_DIRICHLET
        };

    private:
        fluxGridType fgType;
        Grid *grid;
        const len_t STENCIL_WIDTH = 2;
        len_t nnzPerRow = 8*STENCIL_WIDTH-1;
        len_t nr;
        len_t *n1 = nullptr;
        len_t *n2 = nullptr;
        real_t ***deltas = nullptr;
        real_t ***deltas_jac = nullptr;
        real_t *delta_prev;
        real_t *delta_tmp;
        len_t id_unknown;

        // Helper variables that are used in setting coefficients
        // (essentially used like global variables within this class)
        int_t shiftU1, shiftU2, shiftD1;
        real_t xf, x_0, xN;

        adv_bc bc_lower;
        adv_bc bc_upper;
        OptionConstants::adv_jacobian_mode jac_mode = OptionConstants::AD_INTERP_JACOBIAN_LINEAR;

        bool hasBeenInitialized = false;
        bool isFirstRebuild = true;
        bool hasNonTrivialJacobian = false;

        void Deallocate();

        // Returns the "active" index of this interpolation term
        len_t GetIndex(len_t ir, len_t i, len_t j, int_t *N)
        {
            if(fgType == FVM::FLUXGRIDTYPE_RADIAL){
                *N = nr-1;
                return ir;
            } else if(fgType == FVM::FLUXGRIDTYPE_P1){
                *N = n1[ir]-1;
                return i;
            } else {
                *N = n2[ir]-1;
                return j;
            }
        }

        void SetFirstOrderCoefficient(int_t, int_t, const real_t*, real_t, real_t*&, real_t scaleFactor=1.0);
        void SetSecondOrderCoefficient(int_t, int_t, const real_t*, real_t, real_t*&);
        void SetFluxLimitedCoefficient(int_t, int_t, const real_t*, real_t, real_t*&, real_t r=0, real_t psiPrime=0);
        void SetLinearFluxLimitedCoefficient(int_t, int_t, const real_t*, real_t, real_t, real_t*&);
        void SetGPLKScheme(int_t ind, int_t N, const real_t *x, real_t r, real_t alpha, real_t kappa, real_t M, real_t damping, real_t *&deltas);

        std::function<real_t(int_t)> GetYFunc(len_t ir, len_t i, len_t j, FVM::UnknownQuantityHandler *unknowns);
        real_t GetXi(const real_t *x, int_t i, int_t N);
        real_t GetYi(int_t i, int_t N, std::function<real_t(int_t)> y);

        real_t GetFluxLimiterR(int_t ind, int_t N, std::function<real_t(int_t)> y, const real_t *x);

        void SetNNZ(adv_interpolation);

        bool IsFluxLimiterMethod(adv_interpolation method){
            return method==AD_INTERP_TCDF || method==AD_INTERP_OSPRE || method==AD_INTERP_SMART || method==AD_INTERP_MUSCL; 
        }
        bool IsSmoothFluxLimiter(adv_interpolation method) {
            return method==AD_INTERP_TCDF || method==AD_INTERP_OSPRE;
        }
    public:
        AdvectionInterpolationCoefficient(Grid*, fluxGridType,
            adv_bc bc_lower = AD_BC_DIRICHLET,
            adv_bc bc_upper = AD_BC_DIRICHLET
        );
        ~AdvectionInterpolationCoefficient();

        void Allocate();

        void ApplyBoundaryCondition();
        void SetCoefficient(real_t **A, real_t **D=nullptr, UnknownQuantityHandler* unknowns=nullptr, adv_interpolation adv_i=AD_INTERP_CENTRED, real_t damping_factor=1.0);
        bool GridRebuilt();

        void SetUnknownId(len_t id) {id_unknown = id;}
        void SetJacobianMode(OptionConstants::adv_jacobian_mode jac) {jac_mode = jac;}
        void ResetCoefficient();

        void SetBoundaryConditions(adv_bc bc_lower, adv_bc bc_upper){
            this->bc_lower = bc_lower;
            this->bc_upper = bc_upper;
        }
        

        ////////////////////////////////////////////
        // Getters for interpolation coefficients //
        ////////////////////////////////////////////
        real_t ***GetCoefficient(adv_interp_mode interp_mode = AD_INTERP_MODE_FULL) { 
            if(interp_mode == AD_INTERP_MODE_FULL)
                return deltas; 
            else if (interp_mode == AD_INTERP_MODE_JACOBIAN)
                return deltas_jac;
            else
                throw FVMException("Invalid advection interpolation mode requested.");
        }
        real_t *GetCoefficient(len_t ir, len_t i, len_t j,adv_interp_mode interp_mode = AD_INTERP_MODE_FULL) { 
            if(interp_mode == AD_INTERP_MODE_FULL)
                return deltas[ir][j*n1[ir]+i]; 
            else if (interp_mode == AD_INTERP_MODE_JACOBIAN)
                return deltas_jac[ir][j*n1[ir]+i];
            else
                throw FVMException("Invalid advection interpolation mode requested.");
        }
        const real_t GetCoefficient(len_t ir, len_t i, len_t j, len_t n,adv_interp_mode interp_mode = AD_INTERP_MODE_FULL) const { 
            if(interp_mode == AD_INTERP_MODE_FULL)
                return deltas[ir][j*n1[ir]+i][n]; 
            else if (interp_mode == AD_INTERP_MODE_JACOBIAN)
                return deltas_jac[ir][j*n1[ir]+i][n];
            else
                throw FVMException("Invalid advection interpolation mode requested.");
        }

        len_t GetKmin(len_t ind, len_t *n);
        len_t GetKmax(len_t ind, len_t N);
        
        // Returns the number of non-zeroes per row of an advection sterm
        // using this inteprolation coefficient
        len_t GetNNZPerRow(){return nnzPerRow;}
    };
}

#endif/*_DREAM_FVM_ADVECTION_INTERPOLATION_COEFFICIENT_TERM_HPP*/
