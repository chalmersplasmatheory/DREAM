#ifndef _DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP
#define _DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


namespace DREAM::FVM {
    class MomentQuantity : public EquationTerm {
    public:

        // Settings for specifying thresholds to the MomentQuantity.
        enum pThresholdMode {
            P_THRESHOLD_MODE_MIN_MC,
            P_THRESHOLD_MODE_MAX_MC,
            P_THRESHOLD_MODE_MIN_THERMAL,
            P_THRESHOLD_MODE_MAX_THERMAL
        };
    protected:
        real_t *integrand;
        len_t nIntegrand = 0;

        FVM::Grid *fGrid;
        len_t momentId, fId;
        UnknownQuantityHandler *unknowns;

        len_t nnz_per_row;
        void ResetIntegrand(){
            for(len_t i=0; i<this->nIntegrand; i++)
                this->integrand[i] = 0; 
        }
    private:
        real_t pThreshold;
        pThresholdMode pMode;
        bool hasThreshold;

        bool SatisfiesThreshold(len_t ir, len_t i, len_t j);
    public:
        MomentQuantity(
            Grid*, Grid*, len_t, len_t, UnknownQuantityHandler*,
            real_t pThreshold = 0, pThresholdMode pMode = P_THRESHOLD_MODE_MIN_MC
        );
        virtual ~MomentQuantity();

        virtual len_t GetNumberOfNonZerosPerRow() const { return this->nnz_per_row; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const 
            { return GetNumberOfNonZerosPerRow(); }

        virtual bool GridRebuilt() override;

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP*/
