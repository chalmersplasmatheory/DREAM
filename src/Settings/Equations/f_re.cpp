/**
 * Definition of equations relating to the runaway electron
 * distribution function.
 */

#include <string>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/f_re"


/**
 * Define settings for the runaway electron distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_re(Settings *s) {
    DefineDataR(MODULENAME, s, "n0");
    DefineDataR(MODULENAME, s, "T0");
    DefineDataR2P(MODULENAME, s, "init");
}

/**
 * Construct the equation for the runaway electron distribution function.
 * This method is only called if the runaway grid is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equations.
 */
void SimulationGenerator::ConstructEquation_f_re(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    FVM::Operator *eqn     = new FVM::Operator(runawayGrid);

    // Construct kinetic equation
    eqn->AddTerm(new FVM::TransientTerm(
        runawayGrid, eqsys->GetUnknownID(OptionConstants::UQTY_F_RE)
    ));

    eqn->AddTerm(new ElectricFieldTerm(
        runawayGrid, eqsys->GetUnknownHandler(), eqsys->GetRunawayGridType()
    ));
    eqn->AddTerm(new PitchScatterTerm(
        runawayGrid, eqsys->GetRunawayCollisionHandler(), eqsys->GetRunawayGridType(),
        eqsys->GetUnknownHandler()
    ));
    eqn->AddTerm(new EnergyDiffusionTerm(
        runawayGrid, eqsys->GetRunawayCollisionHandler(), eqsys->GetRunawayGridType(),
        eqsys->GetUnknownHandler()
    ));
    eqn->AddTerm(new SlowingDownTerm(
        runawayGrid, eqsys->GetRunawayCollisionHandler(), eqsys->GetRunawayGridType(),
        eqsys->GetUnknownHandler()
    ));

    // The lower p boundary condition depends on whether the hot-tail
    // grid is enabled or not.
    if (eqsys->HasHotTailGrid()) {
        //eqn->A
    // Runaway grid is connected directly to the fluid grid...
    } else {
    }

    // Boundary condition at p=pmax
    eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(runawayGrid, eqn));

    eqsys->SetOperator(OptionConstants::UQTY_F_RE, OptionConstants::UQTY_F_RE, eqn, "3D kinetic equation");
}

