#ifndef _DREAM_SOLVER_NONLINEAR_BACKTRACKER_HPP
#define _DREAM_SOLVER_NONLINEAR_BACKTRACKER_HPP

#include "DREAM/Solver/PhysicalStepAdjuster.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include <petsc.h>
#include <vector>

namespace DREAM {
    class Backtracker : public PhysicalStepAdjuster {
    protected:
        real_t *f0, *f1, *f2;
        real_t *gradf_deltax;
        real_t lambda1, lambda2;

		len_t *monitor, nMonitor=0;

        // Index of non-trivial unknown currently causing
        // 'lambda' to be most limited...
        len_t limitingUnknown=0;

        bool *decreasing;

        // Backtracking iteration counter
        int_t nIteration=0;

        // Number of non-trivial unknowns
        len_t nnu;

        const real_t ALPHA = 1e-4;

        // Methods
        real_t CalculateLambda();
        void EvaluateTargetFunction(Vec&, FVM::BlockMatrix*);
        bool IsDecreasing(bool*);
        void ResetBacktracking();

    public:
        Backtracker(
			std::vector<len_t>&, FVM::UnknownQuantityHandler*,
			IonHandler*, std::vector<std::string>&
		);
        virtual ~Backtracker();

        virtual real_t Adjust(
            len_t, const real_t*, const real_t*,
            Vec&, FVM::BlockMatrix*
        ) override;
    };
}

#endif/*_DREAM_SOLVER_NONLINEAR_BACKTRACKER_HPP*/
