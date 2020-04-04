#ifndef _DREAM_FVM_ADVECTION_TERM_HPP
#define _DREAM_FVM_ADVECTION_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM::FVM {
    class AdvectionTerm : public EquationTerm {
    protected:
        real_t 
        **fr=nullptr, 
        **f1=nullptr, 
        **f2=nullptr;
        bool coefficientsShared = false;

        // Interpolation coefficients
        real_t **deltar=nullptr, **delta1=nullptr, **delta2=nullptr;
        bool interpolationCoefficientsShared = false;

        void AllocateCoefficients();
        void AllocateInterpolationCoefficients();
        void DeallocateCoefficients();
        void DeallocateInterpolationCoefficients();

    public:
        AdvectionTerm(Grid*, bool allocateCoeffs=false);
        ~AdvectionTerm();

        void SetCoefficients(real_t**, real_t**, real_t**);

        // Accessors to advection coefficients
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2)
        { return fr[ir][i2*n1[ir] + i1]; }
        real_t& F1(const len_t ir, const len_t i1, const len_t i2)
        { return f1[ir][i2*(n1[ir]+1) + i1]; }
        real_t& F2(const len_t ir, const len_t i1, const len_t i2)
        { return f2[ir][i2*n1[ir] + i1]; }

        virtual bool GridRebuilt() override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        void SetInterpolationCoefficients(real_t**, real_t**, real_t**);
    };
}

#endif/*_DREAM_FVM_ADVECTION_TERM_HPP*/
