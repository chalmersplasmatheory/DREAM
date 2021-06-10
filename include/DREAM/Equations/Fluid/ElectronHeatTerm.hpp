#ifndef _DREAM_EQUATION_ELECTRON_HEAT_TERM_HPP
#define _DREAM_EQUATION_ELECTRON_HEAT_TERM_HPP

#include "DREAM/EquationSystem.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "DREAM/Equations/Fluid/IonisationHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIIonizLossTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"

/**
 * Implementation of an equation term which represents the total
 * heat of the cold electrons: W_h = (3/2) * n_cold * T_cold
 */
namespace DREAM {
    class ElectronHeatTerm : public FVM::DiagonalQuadraticTerm {
    private:
        real_t constant;
    public:
        ElectronHeatTerm(FVM::Grid* g, const len_t id_other, FVM::UnknownQuantityHandler *u) 
            : FVM::DiagonalQuadraticTerm(g, id_other, u) {
            this->constant = 1.5 * Constants::ec; 
        }
        virtual void SetWeights() override {
            for(len_t i = 0; i<nr; i++)
                weights[i] = constant;
        }
    };
}

#endif/*_DREAM_EQUATION_ELECTRON_HEAT_TERM_HPP*/
