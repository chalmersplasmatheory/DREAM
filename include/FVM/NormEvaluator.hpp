#ifndef _DREAM_FVM_NORM_EVALUATOR_HPP
#define _DREAM_FVM_NORM_EVALUATOR_HPP

#include <vector>
#include "FVM/config.h"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class NormEvaluator {
    protected:
        UnknownQuantityHandler *unknowns;
        const std::vector<len_t> nontrivials;

        // Vector for storing nontrivial-wise norms. Has
        // as many elements as there are non-trivial unknowns
        real_t *retvec;
        len_t nNontrivials;

    public:
        NormEvaluator(UnknownQuantityHandler*, const std::vector<len_t>&);
        virtual ~NormEvaluator();

        const real_t *Norm2(const real_t*);
        void Norm2(const real_t*, real_t*);
        const real_t *Norm2Diff(const real_t*, const real_t*);
        void Norm2Diff(const real_t*, const real_t*, real_t*);
    };
}

#endif/*_DREAM_FVM_NORM_EVALUATOR_HPP*/
