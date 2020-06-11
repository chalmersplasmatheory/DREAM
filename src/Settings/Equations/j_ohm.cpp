/**
 * Definition of equations relating to the radial profile 
 * of ohmic current density j_ohm.
 * The quantity j_ohm corresponds to 
 * j_\Omega / (B/Bmin),
 * which is constant on flux surfaces and proportional to
 * j_ohm ~ \sigma*E_term,
 * where sigma is a neoclassical conductivity
 * containing various geometrical corrections
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentFromConductivityTerm.hpp"

#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

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
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn, "zero");

    // Otherwise, calculate it from conductivity
    } else {
        FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);
        
        eqn2->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()));
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn1, "j_ohm = sigma*Eterm");
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_E_FIELD, eqn2);
        
    }
}

