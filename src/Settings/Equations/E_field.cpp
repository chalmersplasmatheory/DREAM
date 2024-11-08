/**
 * Definition of equations relating to the electric field.
 *
 * Note that we define the electric field as
 *
 *      <E*B>
 *   -----------
 *   sqrt(<B^2>)
 *
 * in DREAM, where 'E' denotes the local electric field, 'B' the
 * local magnetic field and '<X>' denotes the flux-surface average
 * of a quantity 'X'.
 */

#include <string>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Equation/LinearTransientTerm.hpp"
#include "DREAM/Equations/Fluid/AdaptiveHyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/EFieldFromConductivityTerm.hpp"

#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using namespace std;

/**
 * Implementation of a class which represents the Vloop term of the electric
 * field diffusion equation. This operator is to be applied on the electric
 * field, so it rescales the electric field accordingly.
 */
namespace DREAM {
    class VloopTerm : public FVM::DiagonalLinearTerm {
    public:
        VloopTerm(FVM::Grid* g) : FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = 2*M_PI
                    *sqrt(grid->GetRadialGrid()->GetFSA_B2(ir))*grid->GetRadialGrid()->GetBmin(ir);

                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;

                offset += n1[ir]*n2[ir];
            }
        }
    };
}


/**
 * Implementation of a class which represents the dpsi/dt term of the electric
 * field diffusion equation. This term is also multiplied by psi_t' -- the
 * derivative with respect to r of the toroidal flux -- which is how we
 * normalize the equation.
 */
namespace DREAM {
    class dPsiDtTerm : public FVM::LinearTransientTerm {
    public:
        dPsiDtTerm(FVM::Grid* g, const len_t unknownId) : FVM::LinearTransientTerm(g,unknownId){}

        virtual void SetWeights() override {
            len_t offset = 0;
            FVM::RadialGrid *rGrid = grid->GetRadialGrid();
            for (len_t ir = 0; ir < nr; ir++){
                real_t BdotPhi = rGrid->GetBTorG(ir) * rGrid->GetFSA_1OverR2(ir);

                // psit'/VpVol, multiplied by 2*pi
                real_t psitPrimeOverVpVol  = BdotPhi;

                real_t w = -psitPrimeOverVpVol;

                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;

                offset += n1[ir]*n2[ir];
            }
        }
    };
}


#define MODULENAME "eqsys/E_field"
#define MODULENAME_HYPRES "eqsys/psi_p/hyperresistivity"


/**
 * Define options for the electric field module.
 */
void SimulationGenerator::DefineOptions_ElectricField(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED);

    // Prescribed data (in radius+time)
    DefineDataRT(MODULENAME, s, "data");

    // Prescribed initial profile (when evolving E self-consistently)
    DefineDataR(MODULENAME, s, "init");
    

    // Type of boundary condition on the wall
    s->DefineSetting(MODULENAME "/bc/type", "Type of boundary condition to use on the wall for self-consistent E-field", (int_t)OptionConstants::UQTY_V_LOOP_WALL_EQN_SELFCONSISTENT);

    // Minor radius of the wall, defaults to radius of the plasma.
    s->DefineSetting(MODULENAME "/bc/wall_radius", "Minor radius of the inner wall", (real_t) -1);

    // Inverse wall time, defaults to 0 (infinitely conducting wall, 
    // which is equivalent to prescribing V_loop_wall to 0)
    s->DefineSetting(MODULENAME "/bc/inverse_wall_time", "Inverse wall time, representing the conductivity of the first wall", (real_t) 0.0);
    s->DefineSetting(MODULENAME "/bc/R0", "Major radius used to evaluate the external inductance for conductivity of the first wall", (real_t) 0.0);

    // Prescribed data (in time)
    DefineDataT(MODULENAME "/bc", s, "V_loop_wall");

    // Settings for hyperresistive term
    s->DefineSetting(MODULENAME_HYPRES "/mode", "Mode for the hyperresistive term", (int_t)OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_NEGLECT);
	s->DefineSetting(MODULENAME_HYPRES "/grad_j_tot_max", "Maximum current density gradient prior to activation of hyperresistive term", (real_t)0.0);
	s->DefineSetting(MODULENAME_HYPRES "/gradient_normalized", "Flag indicating whether or not 'grad_j_tot_max' is normalized to the average j_tot", (bool)false);
	s->DefineSetting(MODULENAME_HYPRES "/Lambda0", "Value of adaptive diffusion coefficient when enabled", (real_t)0.0);
	s->DefineSetting(MODULENAME_HYPRES "/min_duration", "Minimum duration of the adaptive hyperresistive term", (real_t)0.5e-3);
    DefineDataRT(MODULENAME_HYPRES, s, "Lambda");
}

/**
 * Routines for checking if initial E/j_ohm profiles have been given.
 */
bool SimulationGenerator::HasInitialEfield(EquationSystem *eqsys, Settings *s) {
    real_t *v = LoadDataR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s, "init");
	if (v == nullptr)
		return false;
	else {
		delete [] v;
		return true;
	}
}
bool SimulationGenerator::HasInitialJtot(EquationSystem *eqsys, Settings *s) {
    real_t *v = LoadDataR("eqsys/j_ohm", eqsys->GetFluidGrid()->GetRadialGrid(), s, "init");
	if (v == nullptr)
		return false;
	else {
		delete [] v;
		return true;
	}
}

/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_E_field(
    EquationSystem *eqsys, Settings *s,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    enum OptionConstants::uqty_E_field_eqn type = (enum OptionConstants::uqty_E_field_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED:
            ConstructEquation_E_field_prescribed(eqsys, s);
            break;

        case OptionConstants::UQTY_E_FIELD_EQN_SELFCONSISTENT:
            ConstructEquation_E_field_selfconsistent(eqsys, s, oqty_terms);
            break;

		case OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED_CURRENT:
			ConstructEquation_E_field_prescribed_current(eqsys, s);
			break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_E_FIELD, type
            );
    }
}

/**
 * Construct the equation for a prescribed electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Operator *eqn = new FVM::Operator(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT_intp(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn, "Prescribed");
    // Initial value
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_E_FIELD,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );

    // Set boundary condition psi_wall = 0
    ConstructEquation_psi_wall_zero(eqsys,s);
}

/**
 * Construct the equation for a self-consistent electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_selfconsistent(
    EquationSystem *eqsys, Settings* s,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
	const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);
	const len_t id_j_tot   = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    // Set equations for self-consistent E field evolution
    FVM::Operator *dtTerm = new FVM::Operator(fluidGrid);
    FVM::Operator *Vloop = new FVM::Operator(fluidGrid);

    std::string eqn = "dpsi_p/dt = V_loop";

    // Add transient term -dpsi/dt
    dtTerm->AddTerm(new dPsiDtTerm(fluidGrid, eqsys->GetUnknownID(OptionConstants::UQTY_POL_FLUX)));
    // Add Vloop term
    Vloop->AddTerm(new VloopTerm(fluidGrid));

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

        eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, hypTerm);
        eqn += " + hyperresistivity";
    } else if (hypres_mode == OptionConstants::EQTERM_HYPERRESISTIVITY_MODE_ADAPTIVE) {
		real_t grad_j_tot_max = s->GetReal(MODULENAME_HYPRES "/grad_j_tot_max");
		bool gradient_normalized = s->GetBool(MODULENAME_HYPRES "/gradient_normalized");
		real_t Lambda0 = s->GetReal(MODULENAME_HYPRES "/Lambda0");
		real_t min_duration = s->GetReal(MODULENAME_HYPRES "/min_duration");

		FVM::Operator *hypTerm = new FVM::Operator(fluidGrid);
		AdaptiveHyperresistiveDiffusionTerm *ahrdt = new AdaptiveHyperresistiveDiffusionTerm(
			fluidGrid, eqsys->GetUnknownHandler(),
			grad_j_tot_max, gradient_normalized,
			Lambda0, min_duration
		);

		hypTerm->AddTerm(ahrdt);
		oqty_terms->psi_p_hyperresistive = ahrdt;

		eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, hypTerm);
		eqn += " + hyperresistivity";
	}

    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_POL_FLUX, dtTerm, eqn);
    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, Vloop);
    
	// Specify initialization...
	if (HasInitialEfield(eqsys, s)) {
		/**
		 * Load initial electric field profile.
		 * If the input profile is not explicitly set, then 'SetInitialValue()' is
		 * called with a null-pointer which results in E=0 at t=0
		 */
		real_t *Efield_init = LoadDataR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s, "init");
		eqsys->SetInitialValue(OptionConstants::UQTY_E_FIELD, Efield_init);
		delete [] Efield_init;
	} else if (HasInitialJtot(eqsys, s)) {
		RunawayFluid *REFluid = eqsys->GetREFluid();

		std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_EfieldFromJtot =
			[id_j_tot,REFluid,fluidGrid](FVM::UnknownQuantityHandler *u, real_t *Efield_init) {

			const real_t *j_tot = u->GetUnknownData(id_j_tot);
			const len_t nr = fluidGrid->GetNCells();
			for (len_t ir = 0; ir < nr; ir++) {
				real_t s = REFluid->GetElectricConductivity(ir);
				real_t B = sqrt(fluidGrid->GetRadialGrid()->GetFSA_B2(ir));

				Efield_init[ir] = j_tot[ir]*B / s;
			}
		};

		// Initialize electric field to dummy value to allow RunawayFluid
		// to be initialized first...
		eqsys->SetInitialValue(id_E_field, nullptr);

		eqsys->initializer->AddRule(
			id_E_field,
			EqsysInitializer::INITRULE_EVAL_FUNCTION,
			initfunc_EfieldFromJtot,
			// Dependencies..
			id_j_tot,
			EqsysInitializer::RUNAWAY_FLUID
		);
	} else {
		const string& fromfile = s->GetString("init/fromfile", false);

		if (fromfile == "")
			throw SettingsException(
				"E_field: Self-consistent electric field evolution requested, but neither initial E or j_tot profiles have been specified."
			);
	}

    // Set equation for self-consistent boundary condition
    ConstructEquation_psi_wall_selfconsistent(eqsys,s);
}

/**
 * Construct the equation for the electric field when the ohmic
 * current is prescribed.
 */
void SimulationGenerator::ConstructEquation_E_field_prescribed_current(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Operator *eqnE = new FVM::Operator(eqsys->GetFluidGrid());
    FVM::Operator *eqnj = new FVM::Operator(eqsys->GetFluidGrid());
	FVM::Operator *eqnjre = new FVM::Operator(eqsys->GetFluidGrid());

	eqnE->AddTerm(new FVM::IdentityTerm(eqsys->GetFluidGrid(), -1.0));
	eqnj->AddTerm(
		new EFieldFromConductivityTerm(
			eqsys->GetFluidGrid(), eqsys->GetUnknownHandler(),
			eqsys->GetREFluid()
		)
	);
	eqnjre->AddTerm(
		new EFieldFromConductivityTerm(
			eqsys->GetFluidGrid(), eqsys->GetUnknownHandler(),
			eqsys->GetREFluid(), -1.0
		)
	);

	const len_t id_j_tot   = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
	const len_t id_j_re    = eqsys->GetUnknownID(OptionConstants::UQTY_J_RE);
	const len_t id_E_field = eqsys->GetUnknownID(OptionConstants::UQTY_E_FIELD);

	eqsys->SetOperator(
		id_E_field, id_E_field, eqnE
	);
    eqsys->SetOperator(
		id_E_field, id_j_tot,
		eqnj, "E = (j_tot-j_re) / sigma"
	);
	eqsys->SetOperator(
		id_E_field, id_j_re, eqnjre
	);

	// Initialize electric field to dummy value to allow RunawayFluid
	// to be initialized first...
	eqsys->SetInitialValue(id_E_field, nullptr);

    // Initial value
    eqsys->initializer->AddRule(
        id_E_field,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
		nullptr,
		// Dependencies..
		id_j_tot,
		id_j_re,
		EqsysInitializer::RUNAWAY_FLUID
    );

    // Set boundary condition psi_wall = 0
    ConstructEquation_psi_wall_zero(eqsys,s);
}

