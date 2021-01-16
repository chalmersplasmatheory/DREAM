#ifndef _DREAM_FVM_EQUATION_TERM_HPP
#define _DREAM_FVM_EQUATION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
    class EquationTerm {
    private:
        std::vector<len_t> derivIdsJacobian;
        std::vector<len_t> derivNMultiplesJacobian;

    protected:
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
        // Returns true if derivId is in the derivIdsJacobian vector
        bool HasJacobianContribution(len_t derivId, len_t *nMultiples=nullptr){
            bool hasDerivIdContribution = false;
            for(len_t i_deriv = 0; i_deriv < derivIdsJacobian.size(); i_deriv++){
                if (derivId == derivIdsJacobian[i_deriv]){
                    if(nMultiples != nullptr)
                        *nMultiples = derivNMultiplesJacobian[i_deriv];
                    hasDerivIdContribution = true;
                }
            }
            return hasDerivIdContribution;
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

        virtual len_t GetNumberOfNonZerosPerRow() const = 0;
        virtual len_t GetNumberOfNonZerosPerRow_jac() const {
            return GetNumberOfNonZerosPerRow() + GetNumberOfMultiplesJacobian();
        }

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
