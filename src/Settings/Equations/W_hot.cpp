/**
 * Definition of equations relating to W_hot (the kinetic energy density of hot electrons).
 */

#include <sstream>
#include <iomanip>

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/KineticEnergyFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/ElectronHeatTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/T_hot"


/**
 * Define options for the hot electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_hot(Settings *s) {
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_HOT_MOMENT);
    s->DefineSetting(MODULENAME "/recombination", "Whether to include recombination radiation (true) or ionization energy loss (false)", (bool)false);
    s->DefineSetting(MODULENAME "/halo_region_losses", "Whether to include losses through the halo region (true) or not (false)", (bool)false);
    s->DefineSetting(MODULENAME "/include_NBI", "Whether to include NBI heating term in T_hot evolution", (bool)false);
	s->DefineSetting(MODULENAME "/enabled", "Whether or not T_hot should be enabled", (bool)false);

    // Prescribed data (in radius+time)
    DefineDataRT(MODULENAME, s, "data");

    // Prescribed initial profile (when evolving T self-consistently)
    DefineDataR(MODULENAME, s, "init");

    // Transport settings
    DefineOptions_Transport(MODULENAME, s, false);
}


/**
 * Construct the equation for the hot electron temperature.
 */
void SimulationGenerator::ConstructEquation_T_hot(
	EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
	struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
	FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

	// Construct the equation T_hot = (3/2)*W_hot/n_hot
	FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
	FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

	len_t id_T_hot = eqsys->GetUnknownID(OptionConstants::UQTY_T_HOT);
	len_t id_W_hot = eqsys->GetUnknownID(OptionConstants::UQTY_W_HOT);
	len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);

	Op1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
	Op2->AddTerm(new ElectronHeatTerm(fluidGrid, id_n_hot, eqsys->GetUnknownHandler()));

	eqsys->SetOperator(id_T_hot, id_W_hot, Op1, "W_hot = (2/3)*n_hot*T_hot");
	eqsys->SetOperator(id_T_hot, id_T_hot, Op2);

	// Set initial value for 'T_hot'
	real_t *T0 = LoadDataR("eqsys/f_hot", fluidGrid->GetRadialGrid(), s, "T0");
	eqsys->SetInitialValue(id_T_hot, T0);
	delete [] T0;

	// Construct W_hot equation
	enum OptionConstants::uqty_T_hot_eqn type =
		(enum OptionConstants::uqty_T_hot_eqn)s->GetInteger(MODULENAME "/type");
	
	switch (type) {
		case OptionConstants::UQTY_T_HOT_MOMENT:
			ConstructEquation_W_hot_moment(eqsys, s);
			break;

		case OptionConstants::UQTY_T_HOT_SELF_CONSISTENT:
			ConstructEquation_W_hot_selfconsistent(eqsys, s, adas, nist, amjuel, oqty_terms);
			break;

		default:
			throw SettingsException(
				"Unrecognized equation type for '%s': %d.",
				OptionConstants::UQTY_T_HOT, type
			);
	}
}


/**
 * Construct the equation for the self-consistent hot electron energy density.
 */
void SimulationGenerator::ConstructEquation_W_hot_selfconsistent(
	EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
	struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
	const len_t id_T_hot = eqsys->GetUnknownID(OptionConstants::UQTY_T_HOT);
	const len_t id_W_hot = eqsys->GetUnknownID(OptionConstants::UQTY_W_HOT);
	const len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
	const len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);

	// We use the same logic as for 'T_cold_selfconsistent' to set
	// up 'W_hot_selfconsistent'. Note though that we assign the equation
	// to 'W_hot', whereas for the cold particles it is assigned to 'T_cold'
	// (and not 'W_cold').
	ConstructEquation_T_cold_selfconsistent(
		MODULENAME, eqsys, s, adas, nist, amjuel, &oqty_terms->T_hot,
		id_W_hot, id_T_hot,
		id_W_hot, id_n_hot, id_j_hot,
		true	// <-- "is for T_hot"
	);

	std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_W_hot =
		[id_T_hot, id_n_hot, eqsys](FVM::UnknownQuantityHandler *u, real_t *W_hot_init) {
		
		FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
		const real_t *T = u->GetUnknownData(id_T_hot);
		const real_t *n = u->GetUnknownData(id_n_hot);
		const len_t nr = fluidGrid->GetNr();

		for (len_t ir = 0; ir < nr; ir++) 
			W_hot_init[ir] = 1.5 * Constants::ec * T[ir]*n[ir];
	};

	// Set initial value for 'W_hot'
	eqsys->initializer->AddRule(
		id_W_hot,
		EqsysInitializer::INITRULE_EVAL_FUNCTION,
		initfunc_W_hot,
		// Dependencies
		id_T_hot, id_n_hot
	);
}


/**
 * Construct the equation for the hot kinetic energy density, 'W_hot',
 * and hot electron temperature, 'T_hot'.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_W_hot_moment(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_W_hot = eqsys->GetUnknownID(OptionConstants::UQTY_W_HOT);
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    // Identity part
    FVM::Operator *Op0 = new FVM::Operator(fluidGrid);
    Op0->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    eqsys->SetOperator(id_W_hot, id_W_hot, Op0);

    std::string desc = "Energy moment of f_hot";

    FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
    real_t pThreshold = 0.0;
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        // With collfreq_mode FULL, n_hot is defined as density above some threshold.
        pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        
        std::ostringstream str;
        str <<std::fixed << std::setprecision(3) << pThreshold;
        switch(pMode){
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p>"+str.str()+"*me*c)";
                break;
            } 
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p>"+str.str()+"*pThermal)";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
                desc = "integral(me*c^2(gamma-1)*f_hot), smooth lower limit at p="+str.str()+"*pThermal";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_MC:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p<"+str.str()+"*me*c)";
                break;
            } 
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p<"+str.str()+"*pThermal)";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
                desc = "integral(me*c^2(gamma-1)*f_hot), smooth upper limit at p="+str.str()+"*pThermal";
                break;
            }    
        }
    }
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new KineticEnergyFromDistributionFunction(
            fluidGrid, hottailGrid, id_W_hot, id_f_hot,eqsys->GetUnknownHandler(),pThreshold, pMode)
        );
    eqsys->SetOperator(id_W_hot, id_f_hot, Op1, desc);

    // Set initialization method
    eqsys->initializer->AddRule(
        id_W_hot,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_f_hot
    );
}

