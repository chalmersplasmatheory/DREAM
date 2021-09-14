#ifndef _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_RE_FLUID_TERM_HPP
#define _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_RE_FLUID_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"

/**
 * Implementation of a class which represents the collisional
 * energy transfer from the runaway fluid to the cold electrons,
 * of the form
 *   dW/dt = e*c*Ec*nRE
 * where
 *   Ec = ncold * 4*pi*lnLambda*r0^2*m_e*c^2/e
 * and lnLambda represents the electron-electron Coulomb logarithm 
 * evaluated at a characteristic runaway momentum of 20mc. 
 */

namespace DREAM {
    class CollisionalEnergyTransferREFluidTerm : public FVM::DiagonalComplexTerm {
    private:
        len_t id_ncold, id_Tcold, id_ni;
        CoulombLogarithm *lnLambdaEE;
        real_t constPreFactor; // = 4*pi*r0^2*c*me*c^2

        // momentum at which the e-e coulomb logarithm is evaluated
        const real_t CHARACTERISTIC_RUNAWAY_MOMENTUM = 20;

        // Set weights for this diagonal term
        virtual void SetWeights() override {
            const real_t *ncold = unknowns->GetUnknownData(id_ncold);
            for(len_t ir=0; ir<nr; ir++){
                real_t lnL = lnLambdaEE->evaluateAtP(ir,CHARACTERISTIC_RUNAWAY_MOMENTUM);
                weights[ir] = constPreFactor * lnL * ncold[ir];
            }
        }

        // Set weights jacobian for this diagonal term
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            if(derivId == id_ncold){
                for(len_t ir=0; ir<nr; ir++){
                    real_t lnL = lnLambdaEE->evaluateAtP(ir,CHARACTERISTIC_RUNAWAY_MOMENTUM);
                    diffWeights[ir] = constPreFactor*lnL;
                }
            } else { // derivatives of the coulomb logarithm
                const real_t *ncold = unknowns->GetUnknownData(id_ncold);
                for(len_t n=0; n<nMultiples; n++){
                    for(len_t ir=0; ir<nr; ir++){
                        real_t dLnL = lnLambdaEE->evaluatePartialAtP(ir,CHARACTERISTIC_RUNAWAY_MOMENTUM,derivId,n);
                        diffWeights[nr*n + ir] = constPreFactor * dLnL * ncold[ir];
                    }
                }
            }
        }
    public:
        // Constructor
        CollisionalEnergyTransferREFluidTerm(
            FVM::Grid *g, FVM::UnknownQuantityHandler *u, 
            CoulombLogarithm *lnLEE, real_t scaleFactor=1.0
        ) : FVM::DiagonalComplexTerm(g,u), lnLambdaEE(lnLEE) {
            this->constPreFactor = scaleFactor*4*M_PI*Constants::r0*Constants::r0*Constants::me*Constants::c*Constants::c*Constants::c;
            this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
            this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
            this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

            AddUnknownForJacobian(unknowns, id_ncold);
            AddUnknownForJacobian(unknowns, id_ni);
            AddUnknownForJacobian(unknowns, id_Tcold);
        }

    };
}

#endif /*_DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_RE_FLUID_TERM_HPP*/