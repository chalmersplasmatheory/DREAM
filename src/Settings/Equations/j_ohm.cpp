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
#include "DREAM/Equations/Fluid/AdaptiveHyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/CurrentFromConductivityTerm.hpp"
#include "DREAM/Equations/Fluid/EFieldFromConductivityTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicElectricFieldTerm.hpp"
#include "DREAM/Equations/Fluid/PredictedOhmicCurrentFromDistributionTerm.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

#include <algorithm>

using namespace DREAM;


#define MODULENAME "eqsys/j_ohm"
#define MODULENAME_EFIELD "eqsys/E_field"
#define MODULENAME_HYPRES "eqsys/psi_p/hyperresistivity"


void SimulationGenerator::DefineOptions_j_ohm(Settings *s){
    s->DefineSetting(MODULENAME "/correctedConductivity", "Determines whether to use f_hot's natural ohmic current or the corrected (~Spitzer) value", (int_t) OptionConstants::CORRECTED_CONDUCTIVITY_ENABLED);
    s->DefineSetting(MODULENAME "/conductivityMode", "Determines which formula to use for the conductivity", (int_t) OptionConstants::CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS);

	DefineDataRT(MODULENAME, s, "data");
	DefineDataR(MODULENAME, s, "init");

	s->DefineSetting(MODULENAME "/Ip0", "When specifying an initial j_tot profile, re-scale the profile to get this total plasma current", (real_t)0.0);
	s->DefineSetting(MODULENAME "/j_type", "Current density profile type to initialize from.", (int_t)0.0);
}

/**
 * Construct the equation for the ohmic current density, 'j_ohm',
 * which represents j_\Omega / (B/Bmin) which is constant on the flux surface.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_ohm(
    EquationSystem *eqsys, Settings *s,
	struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    const len_t id_j_ohm = eqsys->GetUnknownID(OptionConstants::UQTY_J_OHM);
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
    const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
        
    bool hottailMode = (eqsys->HasHotTailGrid()) && (eqsys->GetHotTailGrid()->GetNp2(0)==1);
    /** 
     * If using collfreq_mode FULL and not using hot tail mode (Nxi=1), 
     * calculate ohmic current by integrating the distribution, 
     * possibly with conductivity correction
     */
    if(eqsys->HasHotTailGrid() && (collfreq_mode==OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) && !hottailMode){
		FVM::Operator *Op_johm = new FVM::Operator(fluidGrid);
		Op_johm->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
		std::string desc = "j_ohm = ";

        len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
        len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);

        // add total current carried by f_hot
        FVM::Operator *Op_fhot = new FVM::Operator(fluidGrid);
        Op_fhot->AddTerm(new CurrentDensityFromDistributionFunction(
                fluidGrid, eqsys->GetHotTailGrid(), id_j_ohm, id_f_hot,eqsys->GetUnknownHandler()
        ) );
        eqsys->SetOperator(id_j_ohm, id_f_hot, Op_fhot);
        
        // subtract hot current (add with a scaleFactor of -1.0)
        FVM::Operator *Op_jhot = new FVM::Operator(fluidGrid);
        Op_jhot->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        desc += "integral(v_par*f_hot) - j_hot"; 
        eqsys->SetOperator(id_j_ohm, id_j_hot, Op_jhot);
        
        OptionConstants::corrected_conductivity corrCond = (enum OptionConstants::corrected_conductivity)s->GetInteger(MODULENAME "/correctedConductivity");
        if(corrCond == OptionConstants::CORRECTED_CONDUCTIVITY_ENABLED){
            FVM::Operator *Op_E = new FVM::Operator(fluidGrid);
            // add full ohmic current
            Op_E->AddTerm(new CurrentFromConductivityTerm(
                            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler()
            ) );
            // remove predicted numerical ohmic current (add with a scaleFactor of -1.0)
            Op_E->AddTerm(new PredictedOhmicCurrentFromDistributionTerm(
                            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(), eqsys->GetIonHandler(), -1.0
            ) );
            eqsys->SetOperator(id_j_ohm,id_E_field,Op_E);
            desc += " + E*(sigma-sigma_num) [corrected]";

            // Initialization
            eqsys->initializer->AddRule(
                    id_j_ohm,
                    EqsysInitializer::INITRULE_EVAL_EQUATION,
                    nullptr,
                    // Dependencies
                    id_E_field,
                    id_f_hot,
                    EqsysInitializer::RUNAWAY_FLUID
            );

        } else {
            // Initialization
            eqsys->initializer->AddRule(
                    id_j_ohm,
                    EqsysInitializer::INITRULE_EVAL_EQUATION,
                    nullptr,
                    // Dependencies
                    id_f_hot
            );

        }

		eqsys->SetOperator(id_j_ohm,id_j_ohm,Op_johm,desc);    
    // In all other cases, take full spitzer conductivity
    } else { 
        FVM::Operator *Op_E = new FVM::Operator(fluidGrid);
		FVM::Operator *Op_johm = new FVM::Operator(fluidGrid);

		Op_E->AddTerm(new DREAM::OhmicElectricFieldTerm(fluidGrid, -1.0));
        Op_johm->AddTerm(new EFieldFromConductivityTerm(
            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
			1.0, true
        ));

        eqsys->SetOperator(id_j_ohm, id_E_field, Op_E);
		eqsys->SetOperator(id_j_ohm, id_j_ohm, Op_johm);
        std::string desc = "E = j_ohm/sigma"; 

		// Add hyperresistive term
		enum OptionConstants::eqterm_hyperresistivity_mode hypres_mode =
			(enum OptionConstants::eqterm_hyperresistivity_mode)s->GetInteger(MODULENAME_HYPRES "/mode");

		if (hypres_mode == OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_PRESCRIBED) {
			FVM::Interpolator1D *Lambda = LoadDataRT_intp(
				MODULENAME_HYPRES,
				eqsys->GetFluidGrid()->GetRadialGrid(),
				s, "Lambda", true
			);

			FVM::Operator *hypTerm = new FVM::Operator(fluidGrid);
			HyperresistiveDiffusionTerm *hrdt = new HyperresistiveDiffusionTerm(
				fluidGrid, Lambda
			);
			hypTerm->AddTerm(hrdt);
			oqty_terms->psi_p_hyperresistive = hrdt;

			eqsys->SetOperator(id_j_ohm, id_j_tot, hypTerm);
			desc += " + hyperresistivity";
		} else if (
			hypres_mode == OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_ADAPTIVE ||
			hypres_mode == OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL
		) {
			real_t grad_j_tot_max = s->GetReal(MODULENAME_HYPRES "/grad_j_tot_max");
			bool gradient_normalized = s->GetBool(MODULENAME_HYPRES "/gradient_normalized");
			real_t dBB0 = s->GetReal(MODULENAME_HYPRES "/dBB0");
			real_t suppression_level = s->GetReal(MODULENAME_HYPRES "/suppression_level");

			bool localized = (hypres_mode == OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL);

			FVM::Operator *hypTerm = new FVM::Operator(fluidGrid);
			AdaptiveHyperresistiveDiffusionTerm *ahrdt = new AdaptiveHyperresistiveDiffusionTerm(
				fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetIonHandler(),
				grad_j_tot_max, gradient_normalized,
				dBB0, suppression_level, localized
			);

			hypTerm->AddTerm(ahrdt);
			oqty_terms->psi_p_hyperresistive = ahrdt;

			eqsys->SetOperator(id_j_ohm, id_j_tot, hypTerm);
			desc += " + hyperresistivity";
		}

		RunawayFluid *REFluid = eqsys->GetREFluid();
		std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_JOhm =
			[id_j_ohm,id_E_field,fluidGrid,REFluid](FVM::UnknownQuantityHandler *u, real_t *j_ohm_init) {
				const real_t *E_field = u->GetUnknownData(id_j_ohm);
				const len_t nr = fluidGrid->GetNCells();
				for (len_t ir = 0; ir < nr; ir++) {
					real_t s = REFluid->GetElectricConductivity(ir);
					real_t sqrtB2 = std::sqrt(fluidGrid->GetRadialGrid()->GetFSA_B2(ir));

					// j/B = sigma * E / sqrt(<B^2>)
					//   <=>
					// j/(B/Bmin) = sigma * E / sqrt(<B^2>/Bmin^2)
					j_ohm_init[ir] = s*E_field[ir] / sqrtB2;
				}
			};

        // Initialization
        eqsys->initializer->AddRule(
                id_j_ohm,
                EqsysInitializer::INITRULE_EVAL_FUNCTION,
                initfunc_JOhm,
                // Dependencies
                id_E_field,
                EqsysInitializer::RUNAWAY_FLUID
        );
    }
}

