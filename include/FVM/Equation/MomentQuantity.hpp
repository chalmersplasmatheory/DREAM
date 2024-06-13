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
            P_THRESHOLD_MODE_MC=1,
            P_THRESHOLD_MODE_THERMAL=2,
            P_THRESHOLD_MODE_THERMAL_SMOOTH=3
        };

        // Settings for how to integrate over xi
        enum xiIntegralMode {
            XI_MODE_ALL=1,      // Integrate from xi=-1 to xi=1
            XI_MODE_NEG=2,      // Integrate from xi=-1 to xi=0
            XI_MODE_POS=3       // Integrate from xi=0 to xi=+1
        };

		// Direction of p integration (0 to threshold, or threshold to pMax)
		enum pIntegralMode {
			P_INT_MODE_THR2MAX=1,
			P_INT_MODE_ZER2THR=2
		};

    protected:
        real_t *integrand = nullptr;
        real_t *diffIntegrand = nullptr;
        len_t nIntegrand = 0;

        FVM::Grid *fGrid;
        len_t momentId, fId, id_Tcold;
        UnknownQuantityHandler *unknowns;

        //len_t nnz_per_row;

        // Sets all elements of integrand to 0
        void ResetIntegrand(){
            for(len_t i=0; i<this->nIntegrand; i++)
                this->integrand[i] = 0; 
        }
        // Sets all elements of jacobian integrand to 0
        virtual void ResetDiffIntegrand(){
            len_t N = this->nIntegrand*GetMaxNumberOfMultiplesJacobian();
            for(len_t i=0; i<N; i++)
                this->diffIntegrand[i] = 0; 
        }
        virtual void SetDiffIntegrand(len_t){}

        void AllocateDiffIntegrand();

        real_t ThresholdEnvelope(len_t ir, len_t i);
        real_t DiffThresholdEnvelope(len_t ir, len_t i);
		real_t ThresholdValue(len_t ir);
        void AddDiffEnvelope();

		pIntegralMode *pIntMode;

    private:
        real_t pThreshold;
        pThresholdMode pMode;
        xiIntegralMode xiMode;
        bool hasThreshold;

		pIntegralMode defaultPMode = P_INT_MODE_THR2MAX;

        std::vector<len_t> derivIds;
        std::vector<len_t> derivNMultiples;

        static constexpr real_t smoothEnvelopeStepWidth = 2; // for use with P_THRESHOLD_MODE 'SMOOTH'

    public:
        MomentQuantity(
            Grid*, Grid*, len_t, len_t, UnknownQuantityHandler*,
            real_t pThreshold = 0, pThresholdMode pMode = P_THRESHOLD_MODE_MC,
            xiIntegralMode xiMode = XI_MODE_ALL, pIntegralMode defaultPMode=P_INT_MODE_THR2MAX
        );
        virtual ~MomentQuantity();

        // Figure out the maximum number of non-zeros needed per
        // matrix row...
        virtual len_t GetNumberOfNonZerosPerRow() const override {
            len_t nnz_per_row=0;
            for (len_t i = 0; i < fGrid->GetNr(); i++) {
                len_t nc = fGrid->GetMomentumGrid(i)->GetNCells();
                if (nnz_per_row < nc)
                    nnz_per_row = nc;
            }
            return nnz_per_row; 
        }

		static real_t ThresholdValue(real_t pThreshold, pThresholdMode, real_t Tcold);
        static real_t ThresholdEnvelope(len_t i, real_t pThreshold, pThresholdMode, pIntegralMode, MomentumGrid*, real_t T);

        virtual bool GridRebuilt() override;

        virtual bool SetJacobianBlock(const len_t, const len_t, Matrix*, const real_t*) override;
        virtual void SetMatrixElements(Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_FVM_EQUATION_MOMENT_QUANTITY_HPP*/
