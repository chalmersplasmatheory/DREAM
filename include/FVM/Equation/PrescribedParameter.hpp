#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"

#ifndef _DREAM_FVM_EQUATION_PRESCRIBED_PARAMETER_HPP
#define _DREAM_FVM_EQUATION_PRESCRIBED_PARAMETER_HPP

namespace DREAM::FVM {
    class PrescribedParameter : public EquationTerm {
    private:
        real_t *time=nullptr;
        real_t *data=nullptr;

        Interpolator1D *interp=nullptr;

        // Number of time points
        len_t nt;

        real_t currentTime=0;
        const real_t *currentData=nullptr;

        enum Interpolator1D::interp_method interp_method;

    public:
        PrescribedParameter(Grid *g, enum Interpolator1D::interp_method interp=Interpolator1D::INTERP_LINEAR);
        ~PrescribedParameter();

        void DeallocateData();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        const real_t *GetData() { return currentData; }
        void SetData(const len_t, real_t*, real_t*, bool copy=true);
        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_PRESCRIBED_PARAMETER_HPP*/
