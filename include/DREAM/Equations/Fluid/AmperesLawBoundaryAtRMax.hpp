#ifndef _DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_BOUNDARY_AT_RMAX_HPP
#define _DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_BOUNDARY_AT_RMAX_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class AmperesLawBoundaryAtRMax : public BoundaryCondition {
    private:
        const Operator *equation;
        real_t coefficient;
        FVM::Grid *targetGrid;

        /** 
         * sets sign in derivative:
         * +1 for psi(a) and -1 for psi(rmax) 
         */
        real_t scaleFactor;
    public:
        AmperesLawBoundaryAtRMax(DREAM::FVM::Grid*, DREAM::FVM::Grid*, 
            const DREAM::FVM::Operator*, real_t scaleFactor = 1.0);
        virtual ~AmperesLawBoundaryAtRMax();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }

        virtual bool Rebuild(const real_t, UnknownQuantityHandler*) override;

        virtual void AddToJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override;
        virtual void AddToMatrixElements(DREAM::FVM::Matrix*, real_t*) override;
        virtual void AddToVectorElements(real_t*, const real_t*) override;

        // Not implemented (not used)
        virtual void SetJacobianBlock(const len_t, const len_t, DREAM::FVM::Matrix*, const real_t*) override {}
        virtual void SetMatrixElements(DREAM::FVM::Matrix*, real_t*) override {}
        virtual void SetVectorElements(real_t*, const real_t*) override {}
    };
}

#endif/*_DREAM_EQUATIONS_POLOIDAL_FLUX_AMPERES_LAW_BOUNDARY_AT_RMAX_HPP*/
