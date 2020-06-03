#ifndef _DREAM_FVM_ADVECTION_TERM_HPP
#define _DREAM_FVM_ADVECTION_TERM_HPP

namespace DREAM::FVM { class AdvectionTerm; }

#include <softlib/SFile.h>
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

        const real_t *const* GetAdvectionCoeffR() const { return this->fr; }
        const real_t *GetAdvectionCoeffR(const len_t i) const { return this->fr[i]; }
        const real_t *const* GetAdvectionCoeff1() const { return this->f1; }
        const real_t *GetAdvectionCoeff1(const len_t i) const { return this->f1[i]; }
        const real_t *const* GetAdvectionCoeff2() const { return this->f2; }
        const real_t *GetAdvectionCoeff2(const len_t i) const { return this->f2[i]; }

        const real_t *const* GetInterpolationCoeffR() const { return this->deltar; }
        const real_t *GetInterpolationCoeffR(const len_t i) const { return this->deltar[i]; }
        const real_t *const* GetInterpolationCoeff1() const { return this->delta1; }
        const real_t *GetInterpolationCoeff1(const len_t i) const { return this->delta1[i]; }
        const real_t *const* GetInterpolationCoeff2() const { return this->delta2; }
        const real_t *GetInterpolationCoeff2(const len_t i) const { return this->delta2[i]; }

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 7; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void ResetCoefficients();
        void SetCoefficients(real_t**, real_t**, real_t**);

        // Accessors to advection coefficients
        real_t& Fr(const len_t ir, const len_t i1, const len_t i2) {
            if (ir == nr)
                return fr[ir][i2*n1[ir-1] + i1];
            else
                return fr[ir][i2*n1[ir] + i1];
        }
        real_t& F1(const len_t ir, const len_t i1, const len_t i2)
        { return f1[ir][i2*(n1[ir]+1) + i1]; }
        real_t& F2(const len_t ir, const len_t i1, const len_t i2)
        { return f2[ir][i2*n1[ir] + i1]; }

        virtual bool GridRebuilt() override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        void SetInterpolationCoefficients(real_t**, real_t**, real_t**);

        virtual void SaveCoefficientsSFile(const std::string&);
        virtual void SaveCoefficientsSFile(SFile*);
    };
}

#endif/*_DREAM_FVM_ADVECTION_TERM_HPP*/
