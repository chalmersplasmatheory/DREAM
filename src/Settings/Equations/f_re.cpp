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
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticUpper.hpp"
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
	s->DefineSetting(MODULENAME "/boundarycondition", "Type of boundary condition to use at p=pmax.", (int_t)FVM::BC::PXiExternalLoss::BC_PHI_CONST);

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
    len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);

    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    FVM::Operator *eqn     = new FVM::Operator(runawayGrid);

    // Construct kinetic equation
    eqn->AddTerm(new FVM::TransientTerm(
        runawayGrid, id_f_re
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
		len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			runawayGrid, eqsys->GetHotTailGrid(), runawayGrid,
			eqn, id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
		));
        /*eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticUpper(
            eqsys->GetHotTailGrid(), runawayGrid, eqn, id_f_hot, id_f_re
        ));*/
    // Runaway grid is connected directly to the fluid grid...
    } else {
		// TODO add source term for when we're running without a hot-tail grid.
		throw SettingsException(
			"The runaway grid can currently only be run together with a "
			"hot electron grid."
		);
        // Add mirroring boundary condition to p=0...
        //eqn->SetAdvectionBoundaryConditions(FVM::FLUXGRIDTYPE_P1, FVM::AdvectionInterpolationCoefficient::AD_BC_MIRRORED, FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);
    }

    // Boundary condition at p=pmax
	enum FVM::BC::PXiExternalLoss::bc_type bc =
		(enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger(MODULENAME "/boundarycondition");
    eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
		runawayGrid, eqn, id_f_re, id_f_re, nullptr,
		FVM::BC::PXiExternalLoss::BOUNDARY_KINETIC, bc
	));

    eqsys->SetOperator(OptionConstants::UQTY_F_RE, OptionConstants::UQTY_F_RE, eqn, "3D kinetic equation");




    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC){
        // Add avalanche source
        if(eqsys->GetRunawayGridType() != OptionConstants::MOMENTUMGRID_TYPE_PXI)
            throw FVM::FVMException("f_re: Kinetic avalanche source only implemented for p-xi grid.");
        real_t pCutoff = (real_t)s->GetReal("eqsys/n_re/pCutAvalanche");
        FVM::Operator *Op_ava = new FVM::Operator(runawayGrid);
        Op_ava->AddTerm(new AvalancheSourceRP(runawayGrid, eqsys->GetUnknownHandler(), pCutoff, -1.0 ));
        len_t id_n_re = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_N_RE);
        eqsys->SetOperator(id_f_re, id_n_re, Op_ava);
    }




	// Set initial value
	/*if (eqsys->HasHotTailGrid())
		eqsys->SetInitialValue(id_f_re, nullptr);
	else {}
		// TODO
    */
    real_t *n0 = LoadDataR("eqsys/f_hot", runawayGrid->GetRadialGrid(), s, "n0");
    real_t *T0 = LoadDataR("eqsys/f_hot", runawayGrid->GetRadialGrid(), s, "T0");

    ConstructEquation_f_maxwellian(
        OptionConstants::UQTY_F_RE, eqsys, runawayGrid, n0, T0
    );
}

