#ifndef _DREAM_FVM_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_TERM_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class EquationTerm {
    private:
        std::vector<len_t> derivIdsJacobian;
        std::vector<len_t> derivNMultiplesJacobian;

    protected:
        std::string name = "<NOT SET>";

        len_t nr, *n1=nullptr, *n2=nullptr;
        Grid *grid;

        void AllocateMemory();
        void DeallocateMemory();

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this advection term
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, len_t derivId){
            derivIdsJacobian.push_back(derivId);
            derivNMultiplesJacobian.push_back(u->GetUnknown(derivId)->NumberOfMultiples());
        }
        void AddUnknownForJacobian(FVM::UnknownQuantityHandler *u, const char *UQTY){
            AddUnknownForJacobian(u, u->GetUnknownID(UQTY));
        }
        // Returns total number of multiples of the jacobian contributions
        len_t GetNumberOfMultiplesJacobian() const {
            len_t nnz = 0; 
            for(len_t i = 0; i<derivIdsJacobian.size(); i++)
                nnz += derivNMultiplesJacobian[i];
            return nnz;
        }
        len_t GetMaxNumberOfMultiplesJacobian() const {
            len_t maxN = 0; 
            for(len_t i = 0; i<derivIdsJacobian.size(); i++)
                if(derivNMultiplesJacobian[i]>maxN)
                    maxN = derivNMultiplesJacobian[i];
            return maxN;
        }
        

    public:
        EquationTerm(Grid*);
        virtual ~EquationTerm();

        virtual bool GridRebuilt();

        const std::string& GetName() const { return this->name; }
        void SetName(const std::string& n) { this->name = n; }

        virtual len_t GetNumberOfNonZerosPerRow() const = 0;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const {
            return GetNumberOfNonZerosPerRow() + GetNumberOfMultiplesJacobian();
        }

        bool HasJacobianContribution(len_t derivId, len_t *nMultiples=nullptr);

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) = 0;
        /**
         * Sets the block specified by 'uqtyId' and 'derivId' in the
         * given Jacobian matrix. Note that 'uqtyId' and 'derivId' do
         * NOT necessarily correspond to the indices of the matrix block,
         * but should rather be used to identify which unknown parameters
         * should be differentiated, and which should be differentiated
         * _with respect to_.
         *
         * The RETURN value indicates whether or not any elements were
         * set in the jacobian. If 'true', one or more elements were set
         * in the matrix, while 'false' indicates that no elements were
         * set.
         */
        virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, Matrix*, const real_t*) = 0;
        /**
         * Sets the 'matrix' elements which is used when running in 
         * semi-implicit mode. In general, equation terms are non-linear
         * and no natural matrix representation exists; therefore the 
         * default is a fully explicit equation term (by setting the 'rhs'
         * of the equation to the evaluated equation term)
         */
        virtual void SetMatrixElements(Matrix */*mat*/, real_t *rhs) {
            if(rhs != nullptr)
                SetVectorElements(rhs, nullptr);
        }
        virtual void SetVectorElements(real_t*, const real_t*) = 0;

        std::vector<len_t> GetDerivIdsJacobian(){return derivIdsJacobian;}
        std::vector<len_t> GetNMultiplesJacobian(){return derivNMultiplesJacobian;}
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
