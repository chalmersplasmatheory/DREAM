/**
 * Definition of equations relating to the cold electron temperature.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/BindingEnergyTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


#define MODULENAME "eqsys/T_cold"

/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_T_cold(
    EquationSystem *eqsys, Settings *s, ADAS *adas
) {
    enum OptionConstants::uqty_T_cold_eqn type = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        case OptionConstants::UQTY_T_COLD_SELF_CONSISTENT:
            ConstructEquation_T_cold_selfconsistent(eqsys,s,adas);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_T_COLD, type
            );
    }
}

/**
 * Construct the equation for a prescribed temperature.
 */
void SimulationGenerator::ConstructEquation_T_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    eqn->AddTerm(new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp));

    eqsys->SetEquation(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_T_COLD, eqn, "Prescribed");
}


/**
 * Construct the equation for a self-consistent temperature evolution.
 */
void SimulationGenerator::ConstructEquation_T_cold_selfconsistent(
    EquationSystem *eqsys, Settings *s, ADAS *adas
) {
    
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    /**
     * The self-consistent temperature evolution uses an equation
     * for the total cold electron energy W_c (potential + heat) 
     */
    eqsys->SetUnknown(OptionConstants::UQTY_W_COLD, fluidGrid);

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    FVM::Equation *eqn1 = new FVM::Equation(eqsys->GetFluidGrid());
    FVM::Equation *eqn2 = new FVM::Equation(eqsys->GetFluidGrid());
    FVM::Equation *eqn3 = new FVM::Equation(eqsys->GetFluidGrid());

    eqn1->AddTerm(new FVM::TransientTerm(fluidGrid,unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD)) );
    eqn2->AddTerm(new OhmicHeatingTerm(fluidGrid,unknowns));
    eqn3->AddTerm(new RadiatedPowerTerm(fluidGrid,unknowns,eqsys->GetIonHandler(),adas));

    eqsys->SetEquation(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD,eqn1,"dW/dt = j*E - sum(n_e*n_i*L_i) + ...");
    eqsys->SetEquation(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_E_FIELD,eqn2);
    eqsys->SetEquation(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_N_COLD,eqn3);

    ConstructEquation_W_cold(eqsys, s, adas);

}




/**
 * Implementation of a class which represents the Vloop term of the electric field diffusion equation.
 */
namespace DREAM {
    class ElectronHeatTerm : public FVM::DiagonalQuadraticTerm {
    public:
        ElectronHeatTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_COLD),u){}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = 1.5;
        }
    };
}



/**
 * Construct the equation for electron energy content:
 * W_cold = 3n_cold*T_cold/2 + W_binding,
 * where W_binding is the total binding energy of all
 * ions (i.e. the minimum energy required to fully ionise
 * the entire plasma). 
 */
void SimulationGenerator::ConstructEquation_W_cold(
    EquationSystem *eqsys, Settings* /* s */, ADAS */*adas*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    
    FVM::Equation *eqn1 = new FVM::Equation(eqsys->GetFluidGrid());
    FVM::Equation *eqn2 = new FVM::Equation(eqsys->GetFluidGrid());
    FVM::Equation *eqn3 = new FVM::Equation(eqsys->GetFluidGrid());

    
    eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1) );
    eqsys->SetEquation(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_W_COLD, eqn1, "W_c = 3nT/2 + W_bind");

    eqn2->AddTerm(new ElectronHeatTerm(fluidGrid,eqsys->GetUnknownHandler()) );
    eqsys->SetEquation(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_T_COLD, eqn2);
    
    /**
     * TODO: construct equation for W_binding = sum_i(n_i I_i) with I_i the (constant) ionisation energies.
     * Derive TotalBindingEnergyTerm from DiagonalLinearTerm (which must be generalized for ion species) and 
     * load ionisation potentials I_i from suitable atomic database.
     */
    eqn3->AddTerm(new BindingEnergyTerm(fluidGrid,eqsys->GetIonHandler()) );
    eqsys->SetEquation(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_ION_SPECIES, eqn3);
}