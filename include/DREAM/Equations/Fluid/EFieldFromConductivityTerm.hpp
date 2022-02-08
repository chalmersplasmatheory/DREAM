#ifndef _DREAM_EQUATION_FLUID_EFIELD_FROM_CONDUCTIVITY_TERM_HPP
#define _DREAM_EQUATION_FLUID_EFIELD_FROM_CONDUCTIVITY_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implements the term j_ohm/sigma, which gives the electric field
 * strength when the ohmic current is prescribed.
 */
namespace DREAM {
    class EFieldFromConductivityTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
    protected:
        // Set weights for the Jacobian block. Uses differentiated conductivity provided by REFluid. 
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++)
                for (len_t ir = 0; ir < nr; ir++){
					real_t s = REFluid->GetElectricConductivity(ir);
                    real_t dw = REFluid->evaluatePartialContributionConductivity(ir,derivId,n)
								* sqrt(grid->GetRadialGrid()->GetFSA_B2(ir))
                                / (s*s);
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                            diffWeights[offset + i] = dw;
                    offset += n1[ir]*n2[ir];
                }
        }

        // Set weights as the conductivity with a geometric factor 
        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir))
					/ REFluid->GetElectricConductivity(ir);

                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }

    public:
        EFieldFromConductivityTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, RunawayFluid *ref) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref)
        {
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES));
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_EFIELD_FROM_CONDUCTIVITY_TERM_HPP*/
