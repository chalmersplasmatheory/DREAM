/**
 * Definition of equations relating to j_tot.
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/n_cold"

/**
 * Construct the equation for the total plasma current density, 'j_tot'.
 *
 * TODO: Proper equation for j_RE (integral of f_RE + e*c*n_RE_external)
 * 
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_j_tot(
    EquationSystem *eqsys, Settings* /* s */
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    FVM::Equation *eqn0 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn3 = new FVM::Equation(fluidGrid);

    eqn0->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn2->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn3->AddTerm(new FVM::IdentityTerm(fluidGrid, Constants::ec * Constants::c));
    
    eqsys->SetEquation(OptionConstants::UQTY_J_TOT, OptionConstants::UQTY_J_TOT, eqn0, "j_tot = j_ohm + j_hot + e*c*n_RE");
    eqsys->SetEquation(OptionConstants::UQTY_J_TOT, OptionConstants::UQTY_J_OHM, eqn1);
    eqsys->SetEquation(OptionConstants::UQTY_J_TOT, OptionConstants::UQTY_J_HOT, eqn2);
    eqsys->SetEquation(OptionConstants::UQTY_J_TOT, OptionConstants::UQTY_N_RE, eqn3);



}


