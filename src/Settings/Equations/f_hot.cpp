/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

#include <iostream>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;

#define MODULENAME "eqsys/f_hot"


/**
 * Define settings for the hot-tail distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_hot(Settings *s) {
    DefineOptions_f_general(s, MODULENAME);
}

/**
 * Construct the equation for the hot electron distribution
 * function. This function is only called if the hot-tail grid
 * is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_f_hot(
    EquationSystem *eqsys, Settings *s
) {
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
    bool addExternalBC = not eqsys->HasRunawayGrid();
    FVM::Operator *eqn = ConstructEquation_f_general(
        s, MODULENAME, eqsys, id_f_hot, hottailGrid, eqsys->GetHotTailGridType(),
        eqsys->GetHotTailCollisionHandler(), addExternalBC
    );

    // Add kinetic-kinetic boundary condition if necessary...
    if (!addExternalBC) {
		len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			hottailGrid, hottailGrid, eqsys->GetRunawayGrid(), eqn,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
		));
    }
}

