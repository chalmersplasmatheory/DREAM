#ifndef _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP

#include "FVM/Grid/Grid.hpp"

namespace DREAM {
    class RunawaySourceTerm : public FVM::EquationTerm {
    private:
    public:
        RunawaySourceTerm(FVM::Grid*);
        ~RunawaySourceTerm();

        virtual bool GridRebuilt() override;

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t
    };
}

#endif/*_DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP*/
