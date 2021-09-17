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

    public:
        NewtonStepAdjuster(
            std::vector<len_t>& nu, FVM::UnknownQuantityHandler *uqh
        ) : nontrivial_unknowns(nu), unknowns(uqh) {}
        virtual ~NewtonStepAdjuster() {}

        virtual real_t Adjust(
            len_t, const real_t*, const real_t*,
            Vec&, FVM::BlockMatrix*
        ) = 0;
    };
}

#endif/*_DREAM_SOLVER_NONLINEAR_NEWTON_STEP_ADJUSTER_HPP*/

