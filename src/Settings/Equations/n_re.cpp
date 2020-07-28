/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/n_re"


void SimulationGenerator::DefineOptions_n_re(
    Settings *s
) {
    s->DefineSetting(MODULENAME "/avalanche", "Enable/disable secondary (avalanche) generation.", (bool)true);
    s->DefineSetting(MODULENAME "/dreicer", "Model to use for Dreicer generation.", (int_t)OptionConstants::EQTERM_DREICER_MODE_NONE);
    s->DefineSetting(MODULENAME "/Eceff", "Model to use for calculation of the effective critical field.", (int_t)OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL);
}

/**
 * Construct the equation for the runaway electron density, 'n_re'.
 * If the runaway grid is disabled, then n_re is calculated from the
 * particle flux escaping the hot-tail grid, plus any other runaway
 * sources that are enabled.
 */
void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);

    // Add the transient term
    FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
    Op_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));

    // Add avalanche growth rate
    if (s->GetBool(MODULENAME "/avalanche"))
        Op_nRE->AddTerm(new AvalancheGrowthTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );

    // Add Dreicer runaway rate
    enum OptionConstants::eqterm_dreicer_mode dm = 
        (enum OptionConstants::eqterm_dreicer_mode)s->GetInteger(MODULENAME "/dreicer");
    switch (dm) {
        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::CONNOR_HASTIE_NOCORR, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::CONNOR_HASTIE, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::NEURAL_NETWORK, -1.0
            ));
            break;

        default: break;     // Don't add Dreicer runaways
    }

    eqsys->SetOperator(id_n_re, id_n_re, Op_nRE);

    // Add flux from hot tail grid
    if (hottailGrid) {
        FVM::Operator *Op_nRE_fHot = new FVM::Operator(fluidGrid);
        len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

		// BOUNDARY CONDITION
        if (eqsys->GetHotTailGridType() != OptionConstants::MOMENTUMGRID_TYPE_PXI)
            throw NotImplementedException(
                "Currently, the 'DensityFromBoundaryFlux' term only supports "
                "p/xi momentum grids. Hence, you should use a p/xi grid for the "
                "hot-tail distribution function."
            );

		// NOTE We assume that the flux appearing in the equation for 'f_hot'
		// only appears in the (f_hot, f_hot) part of the equation, i.e. in
		// the diagonal block.
		const FVM::Operator *Op = eqsys->GetEquation(id_f_hot)->GetEquation(id_f_hot);

		if (eqsys->HasRunawayGrid()) {
			len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
			// Influx from hot-tail grid (with runaway grid at higher p)
			Op_nRE_fHot->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
				fluidGrid, eqsys->GetHotTailGrid(), eqsys->GetRunawayGrid(),
				Op, id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_DENSITY
			));
		} else {
			// Influx from hot-tail grid (with "nothing" at higher p)
			enum FVM::BC::PXiExternalLoss::bc_type bc =
				(enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger("eqsys/f_hot/boundarycondition");
			Op_nRE_fHot->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
				fluidGrid, Op, id_f_hot, id_n_re, hottailGrid,
				FVM::BC::PXiExternalLoss::BOUNDARY_FLUID, bc
			));
		}

        eqsys->SetOperator(id_n_re, id_f_hot, Op_nRE_fHot, "n_re = [flux from f_hot] + n_re*Gamma_ava");
    } else {
        /*FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
        Op_nRE->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqsys->SetOperator(id_n_re,id_n_re,Op_nRE, "zero");
        eqsys->initializer->AddRule(
            id_n_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );*/
    }

    // Initialize to zero
    eqsys->SetInitialValue(id_n_re, nullptr);
}

