#ifndef _DREAM_EQUATION_FLUID_OHMIC_FIELD_FROM_CONDUCTIVITY_TERM_HPP
#define _DREAM_EQUATION_FLUID_OHMIC_FIELD_FROM_CONDUCTIVITY_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implementation of a class which represents the j/sigma contribution to Ohm's law.
 */
namespace DREAM {
    class OhmicFieldFromConductivityTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        IonHandler *ionHandler;
    protected:
        // Set weights for the Jacobian block. Uses differentiated conductivity provided by REFluid. 
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++)
                for (len_t ir = 0; ir < nr; ir++){
					real_t sigma  = REFluid->GetElectricConductivity(ir);
					real_t dsigma = REFluid->evaluatePartialContributionConductivity(ir,derivId,n);
					real_t avB2 = grid->GetRadialGrid()->GetFSA_B2(ir);
					real_t Bmin   = grid->GetRadialGrid()->GetBmin(ir);

					real_t dw = -avB2/Bmin * dsigma/(sigma*sigma);
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                            diffWeights[offset + i] = dw;
                    offset += n1[ir]*n2[ir];
                }
        }

        // Set weights as the conductivity with a geometric factor 
        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
				real_t sigma = REFluid->GetElectricConductivity(ir);
				real_t avB2 = grid->GetRadialGrid()->GetFSA_B2(ir);
				real_t Bmin = grid->GetRadialGrid()->GetBmin(ir);

				// Bmin comes from the fact that
				//   j in DREAM = j/(B/Bmin)
                real_t w = avB2 / (sigma*Bmin);

                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    public:
        OhmicFieldFromConductivityTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, RunawayFluid *ref, IonHandler *ih) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), ionHandler(ih)
        {
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES));
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_OHMIC_FIELD_FROM_CONDUCTIVITY_TERM_HPP*/
