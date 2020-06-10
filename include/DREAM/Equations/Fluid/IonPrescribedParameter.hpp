#ifndef _ION_PRESCRIBED_PARAMETER_HPP
#define _ION_PRESCRIBED_PARAMETER_HPP

#include "DREAM/IonHandler.hpp"
#include "DREAM/IonInterpolator1D.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    class IonPrescribedParameter : public FVM::EvaluableEquationTerm {
    protected:
        IonHandler *ions;
        // Index of ion species to which this equation term should
        // be applied
        len_t nIons;
        const len_t *ionIndices;
        len_t *Z;
        
        IonInterpolator1D *iondata;

        real_t **currentData;
        real_t lasttime = std::numeric_limits<real_t>::quiet_NaN();

    public:
        IonPrescribedParameter(
            FVM::Grid*, IonHandler*, const len_t, const len_t*, IonInterpolator1D*
        );
        virtual ~IonPrescribedParameter();

        void AllocateData();
        void DeallocateData();

        real_t Evaluate(real_t*, const real_t*, const len_t, const len_t) override;

        virtual bool GridRebuilt() override {
            throw NotImplementedException(
                "'GridRebuilt()' has not been implemented for "
                "'IonPrescribedParameter' or the 'IonInterpolator1D' objects."
            );
        }
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_ION_PRESCRIBED_PARAMETER_HPP*/
