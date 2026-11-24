#ifndef _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_RE_FLUID_TERM_HPP
#define _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_RE_FLUID_TERM_HPP

#include <cmath>
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
 * evaluated at a characteristic runaway momentum of 20mc. For an exterior
 * population, use the outer momentum face and the finite-speed factor c/v
 * of the superthermal free-electron drag used in the kinetic heating moment.
 */

namespace DREAM {
    class CollisionalEnergyTransferREFluidTerm : public FVM::DiagonalComplexTerm {
    private:
        len_t id_ncold, id_Tcold, id_ni;
        CoulombLogarithm *lnLambdaEE;
        const FVM::Grid *exteriorGrid;
        CollisionQuantity::collqty_settings exteriorSettings;
        real_t constPreFactor; // = 4*pi*r0^2*c*me*c^2

        // momentum at which the e-e coulomb logarithm is evaluated
        const real_t CHARACTERISTIC_RUNAWAY_MOMENTUM = 20;

        real_t GetMomentum(len_t ir) const {
            if (exteriorGrid == nullptr)
                return CHARACTERISTIC_RUNAWAY_MOMENTUM;
            const auto *mg = exteriorGrid->GetMomentumGrid(ir);
            return mg->GetP1_f(mg->GetNp1());
        }

        real_t GetPrefactor(real_t p) const {
            return exteriorGrid == nullptr ? constPreFactor : constPreFactor * std::hypot(1.0, p) / p;
        }

        real_t GetLogarithm(len_t ir, real_t p) {
            return exteriorGrid == nullptr ? lnLambdaEE->evaluateAtP(ir,p)
                : lnLambdaEE->evaluateAtP(ir,p,&exteriorSettings);
        }

        // Set weights for this diagonal term
        virtual void SetWeights() override {
            const real_t *ncold = unknowns->GetUnknownData(id_ncold);
            for(len_t ir=0; ir<nr; ir++){
                const real_t p = GetMomentum(ir);
                weights[ir] = GetPrefactor(p) * GetLogarithm(ir,p) * ncold[ir];
            }
        }

        // Set weights jacobian for this diagonal term
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            if(derivId == id_ncold){
                for(len_t ir=0; ir<nr; ir++){
                    const real_t p = GetMomentum(ir);
                    diffWeights[ir] = GetPrefactor(p) * GetLogarithm(ir,p);
                }
            } else { // derivatives of the coulomb logarithm
                const real_t *ncold = unknowns->GetUnknownData(id_ncold);
                for(len_t n=0; n<nMultiples; n++){
                    for(len_t ir=0; ir<nr; ir++){
                        const real_t p = GetMomentum(ir);
                        const real_t dLnL = exteriorGrid == nullptr
                            ? lnLambdaEE->evaluatePartialAtP(ir,p,derivId,n)
                            : lnLambdaEE->evaluatePartialAtP(ir,p,derivId,n,&exteriorSettings);
                        diffWeights[nr*n + ir] = GetPrefactor(p) * dLnL * ncold[ir];
                    }
                }
            }
        }
    public:
        // Constructor
        CollisionalEnergyTransferREFluidTerm(
            FVM::Grid *g, FVM::UnknownQuantityHandler *u, 
            CoulombLogarithm *lnLEE, real_t scaleFactor=1.0, const FVM::Grid *exteriorGrid=nullptr
        ) : FVM::DiagonalComplexTerm(g,u), lnLambdaEE(lnLEE), exteriorGrid(exteriorGrid) {
            exteriorSettings.lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
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
