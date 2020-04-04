#ifndef _DREAM_FVM_DIFFUSION_TERM_HPP
#define _DREAM_FVM_DIFFUSION_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class DiffusionTerm : public EquationTerm {
    protected:
        real_t
            **drr=nullptr,
            **d11=nullptr, **d12=nullptr,
            **d21=nullptr, **d22=nullptr;
        bool coefficientsShared = false;

    public:
        DiffusionTerm(Grid*, bool allocCoefficients=false);
        ~DiffusionTerm();

        void AllocateCoefficients();
        void DeallocateCoefficients();
        void SetCoefficients(real_t**, real_t**, real_t**, real_t**, real_t**);

        // Accessors to diffusion coefficients
        real_t& Drr(const len_t ir, const len_t i1, const len_t i2)
        { return drr[ir][i2*n1[ir] + i1]; }
        real_t& D11(const len_t ir, const len_t i1_f, const len_t i2)
        { return d11[ir][i2*(n1[ir]+1) + i1_f]; }
        real_t& D12(const len_t ir, const len_t i1_f, const len_t i2)
        { return d12[ir][i2*(n1[ir]+1) + i1_f]; }
        real_t& D21(const len_t ir, const len_t i1, const len_t i2_f)
        { return d21[ir][i2_f*n1[ir] + i1]; }
        real_t& D22(const len_t ir, const len_t i1, const len_t i2_f)
        { return d22[ir][i2_f*n1[ir] + i1]; }

        virtual bool GridRebuilt() override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}


#endif/*_DREAM_FVM_DIFFUSION_TERM_HPP*/
