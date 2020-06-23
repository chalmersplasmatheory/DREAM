#ifndef _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_CROSS_HPP
#define _DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_CROSS_HPP

#include "FVM/Equation/BoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class PXiExternalCross : public BoundaryCondition {
    public:
        enum condition_type {
            TYPE_LOWER,     // B.C. at p=p0
            TYPE_UPPER,     // B.C. at p=pmax
            TYPE_DENSITY    // Flow of particles into fluid quantity
        };
    private:
        Grid *lowerGrid, *upperGrid;
        const Operator *equation;
        enum condition_type type;

        /**
         * Advection and diffusion coefficients from the lower grid,
         * averaged over xi and evaluated on the upper grid.
         */
        real_t *upperAp, *upperDpp, *upperDpx;

    public:
        PXiExternalCross(
            DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*,
            const DREAM::FVM::Operator*, enum condition_type
        );
        virtual ~PXiExternalCross();

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

#endif/*_DREAM_FVM_EQUATION_BOUNDARY_CONDITION_P_XI_EXTERNAL_CROSS_HPP*/
