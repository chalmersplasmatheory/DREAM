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
            P_THRESHOLD_MODE_MIN_MC=1,
            P_THRESHOLD_MODE_MIN_THERMAL=2,
            P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH=3,
            P_THRESHOLD_MODE_MAX_MC=4,
            P_THRESHOLD_MODE_MAX_THERMAL=5,
            P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH=6
        };
    protected:
        real_t *integrand = nullptr;
        real_t *diffIntegrand = nullptr;
        len_t nIntegrand = 0;

        FVM::Grid *fGrid;
        len_t momentId, fId;
        UnknownQuantityHandler *unknowns;

        //len_t nnz_per_row;

        // Sets all elements of integrand to 0
        void ResetIntegrand(){
            for(len_t i=0; i<this->nIntegrand; i++)
                this->integrand[i] = 0; 
        }
        // Sets all elements of jacobian integrand to 0
        virtual void ResetDiffIntegrand(){
            for(len_t i=0; i<this->nIntegrand*MaxNMultiple(); i++)
                this->diffIntegrand[i] = 0; 
        }
        virtual void SetDiffIntegrand(len_t){}

        // Adds derivId to list of unknown quantities that contributes to Jacobian of this advection term
        void AddUnknownForJacobian(len_t derivId){
            derivIds.push_back(derivId);
            derivNMultiples.push_back(unknowns->GetUnknown(derivId)->NumberOfMultiples());
        }
        void AllocateDiffIntegrand();

        len_t MaxNMultiple(){
            len_t nMultiples = 0;
            for(len_t it=0; it<derivIds.size(); it++)
                if (derivNMultiples[it]>nMultiples)
                    nMultiples = derivNMultiples[it];
            return nMultiples;
        }
    private:
        real_t pThreshold;
        pThresholdMode pMode;
        bool hasThreshold;

        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;
        real_t smoothEnvelopeStepWidth = 0;

        real_t ThresholdEnvelope(len_t ir, len_t i, len_t j);
        real_t DiffThresholdEnvelope(len_t ir, len_t i, len_t j);
        void AddDiffEnvelope();

    public:
        MomentQuantity(
            Grid*, Grid*, len_t, len_t, UnknownQuantityHandler*,
            real_t pThreshold = 0, pThresholdMode pMode = P_THRESHOLD_MODE_MIN_MC
        );
        virtual ~MomentQuantity();

        // Figure out the maximum number of non-zeros needed per
        // matrix row...
        virtual len_t GetNumberOfNonZerosPerRow() const {
            len_t nnz_per_row=0;
            for (len_t i = 0; i < fGrid->GetNr(); i++) {
                len_t nc = fGrid->GetMomentumGrid(i)->GetNCells();
                if (nnz_per_row < nc)
                    nnz_per_row = nc;
            }
            return nnz_per_row; 
            }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const {
            len_t nnz = GetNumberOfNonZerosPerRow(); 
            for(len_t i = 0; i<derivIds.size(); i++)
                nnz += derivNMultiples[i];
            return nnz;
            }

        virtual bool GridRebuilt() override;
        //virtual void Rebuild(const real_t, const real_t, UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP*/
