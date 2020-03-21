#ifndef _TQS_FVM_EQUATION_TERM_HPP
#define _TQS_FVM_EQUATION_TERM_HPP

#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/Matrix.hpp"

namespace TQS::FVM {
    class EquationTerm {
    protected:
        RadialGrid *grid;

        // Interpolation coefficients
        real_t **deltar=nullptr, **delta1=nullptr, **delta2=nullptr;
        bool interpolationCoeffsShared = false;

        void AllocateInterpolationCoefficients();
        void DeallocateInterpolationCoefficients();

    public:
        EquationTerm(RadialGrid*, bool allocInterpolationCoeffs=true);
        ~EquationTerm();

        virtual bool GridRebuilt();

        virtual void Rebuild(const real_t) = 0;
        virtual void SetMatrixElements(Matrix*) = 0;

        void SetInterpolationCoefficients(real_t**, real_t**, real_t**);
    };
}

#endif/*_TQS_FVM_EQUATION_TERM_HPP*/
