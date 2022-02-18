#ifndef _ION_SOURCE_TERM_HPP
#define _ION_SOURCE_TERM_HPP

#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/Fluid/IonPrescribedParameter.hpp"
#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Equation/EvaluableEquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
    class IonSourceTerm : public IonPrescribedParameter {
    public:
        IonSourceTerm(
            FVM::Grid*, IonHandler*, const len_t, const len_t*,
			MultiInterpolator1D*
        );
        virtual ~IonSourceTerm();

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_ION_SOURCE_TERM_HPP*/
