#ifndef _DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP
#define _DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implementation of a class which represents the Gamma_ava*n_re contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
namespace DREAM {
    class AvalancheGrowthTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        len_t id_n_re;
        real_t scaleFactor;
    protected:
        // Set weights for the Jacobian block. Uses differentiated growth rate provided by REFluid. 
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            real_t *dGamma = REFluid->evaluatePartialContributionAvalancheGrowthRate(derivId);

            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++){
                for (len_t ir = 0; ir < nr; ir++){
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                        diffWeights[offset + i] = scaleFactor*dGamma[offset + i];
                    offset += n1[ir]*n2[ir];
                }
            }
            delete [] dGamma;
        }

        // Set weights as the avalanche growth rate
        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t Gamma_ava = REFluid->GetAvalancheGrowthRate(ir);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = scaleFactor*Gamma_ava;
                offset += n1[ir]*n2[ir];
            }
        }
    public:
        AvalancheGrowthTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, real_t scaleFactor = 1.0) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), scaleFactor(scaleFactor)
        {
            id_n_re = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);

            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP*/
