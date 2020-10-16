/**
 * Definition of equations relating to the runaway electron
 * distribution function.
 */

#include <string>
#include "DREAM/DREAMException.hpp"
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


#define MODULENAME "eqsys/f_re"


/**
 * Define settings for the runaway electron distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_re(Settings *s) {
    DefineOptions_f_general(s, MODULENAME);
}

/**
 * Construct the equation for the runaway electron distribution function.
 * This method is only called if the runaway grid is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equations.
 */
void SimulationGenerator::ConstructEquation_f_re(
    EquationSystem *eqsys, Settings *s, struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
    bool addExternalBC = true;
    bool addInternalBC = false;
    FVM::Operator *eqn = ConstructEquation_f_general(
        s, MODULENAME, eqsys, id_f_re, runawayGrid, eqsys->GetRunawayGridType(),
        eqsys->GetRunawayCollisionHandler(), addExternalBC, addInternalBC,
        &oqty_terms->f_re_advective_bc, &oqty_terms->f_re_diffusive_bc
    );

    // Add kinetic-kinetic boundary condition if necessary...
    if (eqsys->HasHotTailGrid()) {
		len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        const FVM::Operator *eqn_f_hot = eqsys->GetOperator(id_f_hot, id_f_hot);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			runawayGrid, eqsys->GetHotTailGrid(), runawayGrid, eqn_f_hot,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
		));
    } else
        throw DREAMException("f_re: Cannot run with RE distribution and without HOT distribution.");
}

