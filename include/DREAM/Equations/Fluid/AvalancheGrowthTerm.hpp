#ifndef _DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP
#define _DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
/**
 * Implementation of a class which represents the Gamma_ava*n_re contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
namespace DREAM {
    class AvalancheGrowthTerm : public FVM::DiagonalComplexTerm, public RunawaySourceTerm {
    private:
        RunawayFluid *REFluid;
        len_t id_n_re;
        real_t scaleFactor;
        real_t *dGamma = nullptr;
        len_t nr_tmp=0;

    protected:
        // Set weights for the Jacobian block. Uses differentiated growth rate provided by REFluid. 
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            AllocateDGamma();
            REFluid->evaluatePartialContributionAvalancheGrowthRate(dGamma, derivId);
            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++){
                for (len_t ir = 0; ir < nr; ir++){
                    const len_t xiIndex = this->GetXiIndexForEDirection(ir);

                    diffWeights[offset + n1[ir]*xiIndex] = scaleFactor*dGamma[ir];
                    offset += n1[ir]*n2[ir];
                }
            }
        }

        // Set weights as the avalanche growth rate
        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                const real_t *GammaAva = REFluid->GetAvalancheGrowthRate();
                const len_t xiIndex = this->GetXiIndexForEDirection(ir);

                weights[offset + n1[ir]*xiIndex] = scaleFactor*GammaAva[ir];
                offset += n1[ir]*n2[ir];
            }
        }

        // if nr has changed, (re)allocate dGamma
        void AllocateDGamma(){
            if(nr_tmp != this->nr){
                if(dGamma != nullptr)
                    delete [] dGamma;
                dGamma = new real_t[this->nr];
                nr_tmp = this->nr;
            }
        }

    public:
        AvalancheGrowthTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, real_t scaleFactor = 1.0) 
            : FVM::DiagonalComplexTerm(g,u), RunawaySourceTerm(g,u), REFluid(ref), scaleFactor(scaleFactor)
        {
            id_n_re = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
        }

        ~AvalancheGrowthTerm(){
            if(dGamma != nullptr)
                delete [] dGamma;
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP*/
