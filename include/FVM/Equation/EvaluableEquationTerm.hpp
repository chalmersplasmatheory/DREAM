#ifndef _DREAM_FVM_EQUATION_EVALUABLE_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_EVALUABLE_EQUATION_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class EvaluableEquationTerm : public EquationTerm {
    public:
        EvaluableEquationTerm(Grid *g) : EquationTerm(g) {}

        virtual void Evaluate(real_t*, const real_t*, const len_t, const len_t) = 0;
    };
}

#endif/*_DREAM_FVM_EQUATION_EVALUABLE_EQUATION_TERM_HPP*/
