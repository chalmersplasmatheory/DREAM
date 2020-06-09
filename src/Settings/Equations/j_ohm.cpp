/**
 * Definition of equations relating to the radial profile 
 * of ohmic current density j_ohm.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include <algorithm>

using namespace DREAM;

#define MODULENAME "eqsys/j_ohm"


/**
 * Construct the equation for the ohmic current density, 'j_ohm',
 * which represents j_\Omega / (B/Bmin) which is constant on the flux surface.
 * This is zero when the hot tail grid uses collfreqmode FULL,
 * in which case the ohmic current is part of f_hot. 
 * The ohmic current j_ohm is calculated from a semi-analytical 
 * electric conductivity formula.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_ohm(
    EquationSystem *eqsys, Settings *s 
) {

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();

    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");


    // If collfreqmode is FULL, the ohmic current is naturally captured in j_hot and we set j_ohm=0.
    if (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid));

        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn, "zero");

    // Otherwise, calculate it from conductivity
    } else {
        RunawayFluid *REFluid = eqsys->GetREFluid();
        IonHandler *ionHandler = eqsys->GetIonHandler();

        FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);
        
        // weightFunc represents sigma * sqrt(<B^2/Bmin^2>)
        std::function<real_t(len_t,len_t,len_t)> weightFunc = [fluidGrid,REFluid,ionHandler](len_t ir,len_t, len_t)
            {return sqrt(fluidGrid->GetRadialGrid()->GetFSA_B2(ir)) * REFluid->evaluateElectricalConductivity(ir, ionHandler->evaluateZeff(ir));};

        eqn2->AddTerm(new FVM::WeightedIdentityTerm(fluidGrid, &weightFunc));
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn1, "Ohmic current from conductivity");
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_E_FIELD, eqn2);
        
    }
}

