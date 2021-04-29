#ifndef _DREAM_FVM_DIAGONAL_TERM_HPP
#define _DREAM_FVM_DIAGONAL_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM::FVM {
    class DiagonalTerm : public EquationTerm {
    private:    
        virtual void DeallocateWeights();

    protected:
        bool hasBeenInitialized = false;
        virtual void AllocateWeights();
        virtual void AllocateDiffWeights(){}
        real_t *weights = nullptr;

        virtual len_t GetNumberOfWeightsElements()
            {return grid->GetNCells();}
        virtual void InitializeWeights();
        virtual bool TermDependsOnUnknowns() = 0; // determines whether weights should be set at every Rebuild or just on GridRebuilt
        virtual void SetWeights() = 0;
        virtual bool AddWeightsJacobian(const len_t, const len_t, Matrix*, const real_t*) = 0;
    public:
        DiagonalTerm(Grid*);
        ~DiagonalTerm();
        
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }

        virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;
        virtual bool GridRebuilt() override;
        virtual bool SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_DIAGONAL_TERM_HPP*/
