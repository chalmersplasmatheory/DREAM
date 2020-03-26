#ifndef _DREAM_FVM_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM::FVM {
    class EquationTerm {
    protected:
        len_t nr, *n1=nullptr, *n2=nullptr;
        Grid *grid;

        // Interpolation coefficients
        real_t **deltar=nullptr, **delta1=nullptr, **delta2=nullptr;
        bool interpolationCoeffsShared = false;

        void AllocateInterpolationCoefficients();
        void DeallocateInterpolationCoefficients();

    public:
        EquationTerm(Grid*, bool allocInterpolationCoeffs=true);
        ~EquationTerm();

        virtual bool GridRebuilt();

        virtual void Rebuild(const real_t) = 0;
        virtual void SetMatrixElements(Matrix*) = 0;

        void SetInterpolationCoefficients(real_t**, real_t**, real_t**);
    };
}

#endif/*_DREAM_FVM_EQUATION_TERM_HPP*/
