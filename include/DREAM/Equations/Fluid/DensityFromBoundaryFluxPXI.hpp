#ifndef _DREAM_DENSITY_FROM_BOUNDARY_FLUX_HPP
#define _DREAM_DENSITY_FROM_BOUNDARY_FLUX_HPP

#include <functional>
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class DensityFromBoundaryFluxPXI : public FVM::EquationTerm {
    private:
        FVM::Grid *distributionGrid;
        const FVM::Operator *equation;

        void __SetElements(std::function<void(const len_t, const len_t, const real_t)>);
    public:
        DensityFromBoundaryFluxPXI(FVM::Grid*, FVM::Grid*, const FVM::Operator*);
        ~DensityFromBoundaryFluxPXI();

        virtual len_t GetNumberOfNonZerosPerRow() const override;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override;

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_DENSITY_FROM_BOUNDARY_FLUX_HPP*/
