#ifndef _DREAM_FVM_BOUNDARY_CONDITION_HPP
#define _DREAM_FVM_BOUNDARY_CONDITION_HPP

#include <string>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM::BC {
    class BoundaryCondition {
    protected:
        Grid *grid;
        std::string name;

    public:
        BoundaryCondition(Grid *g) : grid(g) {};
        virtual ~BoundaryCondition() {}

        const std::string& GetName() const { return this->name; }
        void SetName(const std::string& n) { this->name = n; }

        // By default, we assume that these do not add any
        // additional non-zero elements (which did not already
        // exist after applying the regular operator). This must
        // of course be corrected for in case it is not true.
        virtual len_t GetNumberOfNonZerosPerRow() const { return 0; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return this->GetNumberOfNonZerosPerRow(); }

        virtual bool GridRebuilt() { return false; }

        virtual bool Rebuild(const real_t t, UnknownQuantityHandler*) = 0;

        virtual bool AddToJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) = 0;
        virtual bool SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) = 0;
        virtual void AddToMatrixElements(Matrix*, real_t*) =0;
        virtual void SetMatrixElements(Matrix*, real_t*) = 0;
        virtual void AddToVectorElements(real_t*, const real_t*) = 0;
        virtual void SetVectorElements(real_t*, const real_t*) = 0;
    };
}

#endif/*_DREAM_FVM_BOUNDARY_CONDITION_HPP*/
