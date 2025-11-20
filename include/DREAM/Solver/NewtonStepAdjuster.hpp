#ifndef _DREAM_SOLVER_NONLINEAR_NEWTON_STEP_ADJUSTER_HPP
#define _DREAM_SOLVER_NONLINEAR_NEWTON_STEP_ADJUSTER_HPP

#include <petsc.h>
#include <vector>
#include "FVM/BlockMatrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class NewtonStepAdjuster {
    protected:
        std::vector<len_t> nontrivial_unknowns;
        FVM::UnknownQuantityHandler *unknowns;
		len_t vector_size;

		Vec *residual;
		FVM::BlockMatrix *jacobian;

		real_t *x0=nullptr;
		const real_t *dx=nullptr;

		void UpdateSolution(real_t*, const real_t);

    public:
        NewtonStepAdjuster(
            std::vector<len_t>&, FVM::UnknownQuantityHandler*, const len_t
        );
        virtual ~NewtonStepAdjuster();

		virtual bool AdjustmentNeeded(const len_t, Vec&, FVM::BlockMatrix*) = 0;
		virtual void AdjustSolution(const len_t, real_t*) = 0;
		virtual void Reset(Vec&, FVM::BlockMatrix*) = 0;

		virtual void SetX0(const len_t, const real_t*, const real_t*);

		virtual real_t GetCurrentDamping() = 0;
    };
}

#endif/*_DREAM_SOLVER_NONLINEAR_NEWTON_STEP_ADJUSTER_HPP*/

