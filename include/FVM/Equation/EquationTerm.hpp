#ifndef _DREAM_FVM_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"

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

        virtual void Rebuild(const real_t) = 0;
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
