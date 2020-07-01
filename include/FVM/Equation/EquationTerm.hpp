#ifndef _DREAM_FVM_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class EquationTerm {
    protected:
        len_t nr, *n1=nullptr, *n2=nullptr;
        Grid *grid;

        void AllocateMemory();
        void DeallocateMemory();

    public:
        EquationTerm(Grid*);
        virtual ~EquationTerm();

        virtual bool GridRebuilt();

        virtual len_t GetNumberOfNonZerosPerRow() const = 0;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const = 0;

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) = 0;
        /**
         * Sets the block specified by 'uqtyId' and 'derivId' in the
         * given Jacobian matrix. Note that 'uqtyId' and 'derivId' do
         * NOT necessarily correspond to the indices of the matrix block,
         * but should rather be used to identify which unknown parameters
         * should be differentiated, and which should be differentiated
         * _with respect to_.
         */
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*) = 0;
        virtual void SetMatrixElements(Matrix*, real_t*) = 0;
        virtual void SetVectorElements(real_t*, const real_t*) = 0;
    };

    class EquationTermException : public FVMException {
    public:
        template<typename ... Args>
        EquationTermException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("EquationTerm");
        }
    };
}

#endif/*_DREAM_FVM_EQUATION_TERM_HPP*/
