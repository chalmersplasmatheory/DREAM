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
        real_t *dGamma;
        real_t *gamma;
    protected:
        // Set weights for the Jacobian block. Uses differentiated growth rate provided by REFluid. 
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            //real_t *dGamma = REFluid->evaluatePartialContributionAvalancheGrowthRate(derivId);
            REFluid->evaluatePartialContributionAvalancheGrowthRate(dGamma, derivId);

            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++){
                for (len_t ir = 0; ir < nr; ir++){
                    for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                        diffWeights[offset + i] = scaleFactor*dGamma[offset + i];
                    offset += n1[ir]*n2[ir];
                }
            }
        }

        // Set weights as the avalanche growth rate
        virtual void SetWeights() override {
            len_t offset = 0;
//            real_t *nRE = unknowns->GetUnknownData(id_n_re);
            for (len_t ir = 0; ir < nr; ir++){
                REFluid->SetAvalancheGrowthRate(gamma);
//                real_t sgn_nre = (nRE[ir]>0) - (nRE[ir]<0); 
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = scaleFactor*gamma[offset + i];
                offset += n1[ir]*n2[ir];
            }
        }

	/**
	 * Allocate memory for the runaway rate.
	 */
	void AllocateGamma() {
	    this->gamma = new real_t[this->grid->GetNr()];
	    this->dGamma = new real_t[this->grid->GetNr()];
	}

	/**
	 * Free memory for the runaway rate.
	 */
	void DeallocateGamma() {
	    delete [] this->gamma;
	    delete [] this->dGamma;
	}

	/**
	 * Method called when the grid has been rebuilt.
	 */
	virtual bool GridRebuilt() override {
	    DeallocateGamma();
	    AllocateGamma();
	    return FVM::DiagonalComplexTerm::GridRebuilt();
	}
    public:
        AvalancheGrowthTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, real_t scaleFactor = 1.0) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), scaleFactor(scaleFactor)
        {
            id_n_re = this->unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
            AllocateGamma();
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD));
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
        }
        ~AvalancheGrowthTerm(){
            DeallocateGamma();
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_AVALANCHE_GROWTH_TERM_HPP*/
