/**
 * Definition of equations relating to j_tot, 
 * representing the total plasma current:
 *  j_tot = j_ohm + j_hot + j_RE
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_tot"

/**
 * Define options for the electric field module.
 */
void SimulationGenerator::DefineOptions_j_tot(Settings*){
    //s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED);
}



/**
 * Construct the equation for the total plasma current density, 'j_tot',
 * and for the total plasma current 'I_p'.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_j_tot(
    EquationSystem *eqsys, Settings *s
) {
	enum OptionConstants::uqty_E_field_eqn type =
		(enum OptionConstants::uqty_E_field_eqn)s->GetInteger("eqsys/E_field/type");
	
	if (type == OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED_CURRENT)
		ConstructEquation_j_tot_prescribed(eqsys, s);
	else
		ConstructEquation_j_tot_consistent(eqsys, s);

	ConstructEquation_Ip(eqsys);
}

void SimulationGenerator::ConstructEquation_j_tot_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
	FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
	FVM::RadialGrid *rGrid = fluidGrid->GetRadialGrid();

	FVM::Operator *eqn = new FVM::Operator(fluidGrid);

	struct dream_2d_data *dd = LoadDataRT("eqsys/j_ohm", rGrid, s);

	// Re-scale to get right initial plasma current (Ip)?
	if (s->GetReal("eqsys/j_ohm/Ip0") != 0) {
		const real_t Ip0 = s->GetReal("eqsys/j_ohm/Ip0");
		real_t Ipj = TotalPlasmaCurrentFromJTot::EvaluateIp(rGrid, dd->x);

		const real_t f = Ip0 / Ipj;
		const len_t N = dd->nt * dd->nr;
		for (len_t i = 0; i < N; i++)
			dd->x[i] *= f;
	}

	FVM::Interpolator1D *interp = new FVM::Interpolator1D(dd->nt, dd->nr, dd->t, dd->x, dd->interp);
	FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
	eqn->AddTerm(pp);

	eqsys->SetOperator(id_j_tot, id_j_tot, eqn, "Prescribed");
	// Initial value
	eqsys->initializer->AddRule(
		id_j_tot,
		EqsysInitializer::INITRULE_EVAL_EQUATION
	);
}

void SimulationGenerator::ConstructEquation_j_tot_consistent(
	EquationSystem *eqsys, Settings *s
) {
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
    const len_t id_j_ohm = eqsys->GetUnknownID(OptionConstants::UQTY_J_OHM);
    const len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
    const len_t id_j_re  = eqsys->GetUnknownID(OptionConstants::UQTY_J_RE);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    FVM::Operator *eqn0 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn1 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn2 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn3 = new FVM::Operator(fluidGrid);

    
    eqn0->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn2->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn3->AddTerm(new FVM::IdentityTerm(fluidGrid));
    
    eqsys->SetOperator(id_j_tot, id_j_tot, eqn0, "j_tot = j_ohm + j_hot + j_re");
    eqsys->SetOperator(id_j_tot, id_j_ohm, eqn1);
    eqsys->SetOperator(id_j_tot, id_j_hot, eqn2);
    eqsys->SetOperator(id_j_tot, id_j_re,  eqn3);
    
    // Initialization
	if (HasInitialJtot(eqsys, s)) {
		real_t *jtot_init = LoadDataR("eqsys/j_ohm", eqsys->GetFluidGrid()->GetRadialGrid(), s, "init");

		// Re-scale to get right plasma current (Ip)?
		if (s->GetReal("eqsys/j_ohm/Ip0") != 0) {
			FVM::RadialGrid *rGrid = fluidGrid->GetRadialGrid();
			const len_t nr = rGrid->GetNr();
			const real_t Ip0 = s->GetReal("eqsys/j_ohm/Ip0");
			real_t Ipj = TotalPlasmaCurrentFromJTot::EvaluateIp(rGrid, jtot_init);

			for (len_t ir = 0; ir < nr; ir++)
				jtot_init[ir] *= Ip0 / Ipj;
		}

		eqsys->SetInitialValue(OptionConstants::UQTY_J_TOT, jtot_init);
		delete [] jtot_init;
	} else {
		eqsys->initializer->AddRule(
			id_j_tot,
			EqsysInitializer::INITRULE_EVAL_EQUATION,
			nullptr,
			// Dependencies
			id_j_ohm, id_j_hot, id_j_re
		);
	}
}

void SimulationGenerator::ConstructEquation_Ip(
	EquationSystem *eqsys
) {
    const len_t id_I_p   = eqsys->GetUnknownID(OptionConstants::UQTY_I_P);
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    /**
     * Set equation for the total plasma current I_p
     * (as an integral over j_tot).
     */
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::Operator *eqn_Ip1 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Ip2 = new FVM::Operator(scalarGrid);
    
    eqn_Ip1->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqn_Ip2->AddTerm(
		new TotalPlasmaCurrentFromJTot(
			scalarGrid, fluidGrid,
			eqsys->GetUnknownHandler(), id_j_tot
		)
	);

    eqsys->SetOperator(id_I_p, id_I_p,   eqn_Ip1, "Ip = integral(j_tot)");
    eqsys->SetOperator(id_I_p, id_j_tot, eqn_Ip2);

    eqsys->initializer->AddRule(
        id_I_p,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_j_tot
    );
}


