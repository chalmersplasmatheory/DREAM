/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/Equations/Fluid/KineticEquationTermIntegratedOverMomentum.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Equations/TransportPrescribed.hpp"
#include "DREAM/IO.hpp"
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
    s->DefineSetting(MODULENAME "/avalanche", "Model to use for secondary (avalanche) generation.", (int_t) OptionConstants::EQTERM_AVALANCHE_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/pCutAvalanche", "Minimum momentum to which the avalanche source is applied", (real_t) 0.0);
    s->DefineSetting(MODULENAME "/dreicer", "Model to use for Dreicer generation.", (int_t)OptionConstants::EQTERM_DREICER_MODE_NONE);
    s->DefineSetting(MODULENAME "/Eceff", "Model to use for calculation of the effective critical field.", (int_t)OptionConstants::COLLQTY_ECEFF_MODE_FULL);

    s->DefineSetting(MODULENAME "/adv_interp/r", "Type of interpolation method to use in r-component of advection term of kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(MODULENAME "/adv_interp/r_jac", "Type of interpolation method to use in the jacobian of the r-component of advection term of kinetic equation.", (int_t)OptionConstants::AD_INTERP_JACOBIAN_LINEAR);
    s->DefineSetting(MODULENAME "/adv_interp/fluxlimiterdamping", "Underrelaxation parameter that may be needed to achieve convergence with flux limiter methods", (real_t) 1.0);

    DefineOptions_Transport(MODULENAME, s, false);

    s->DefineSetting(MODULENAME "/compton/mode", "Model to use for Compton seed generation.", (int_t) OptionConstants::EQTERM_COMPTON_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/compton/flux", "Gamma ray photon flux (m^-2 s^-1).", (real_t) 0.0);

    s->DefineSetting(MODULENAME "/tritium", "Indicates whether or not tritium decay RE generation should be included.", (bool)false);

    s->DefineSetting(MODULENAME "/hottail", "Model to use for hottail runaway generation.", (int_t) OptionConstants::EQTERM_HOTTAIL_MODE_DISABLED);

    // Prescribed initial profile
    DefineDataR(MODULENAME, s, "init");

}

/**
 * Construct the equation for the runaway electron density, 'n_re'.
 * It is given by the particle flux escaping the hot-tail grid, 
 * plus any other runaway sources that are enabled.
 */
void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings *s,
    struct OtherQuantityHandler::eqn_terms *oqty_terms,
    FVM::Operator *transport_fre
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();

    len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    len_t id_n_tot = eqsys->GetUnknownID(OptionConstants::UQTY_N_TOT);
    len_t id_n_i   = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    // Add the transient term
    FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_n_tot = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_n_i = new FVM::Operator(fluidGrid);
    Op_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));

    std::string desc_sources = ""; 

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
		const FVM::Operator *Op = eqsys->GetEquation(id_f_hot)->GetOperator(id_f_hot);

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
				fluidGrid, Op, id_f_hot, hottailGrid,
				FVM::BC::PXiExternalLoss::BOUNDARY_FLUID, bc
			));
		}
        desc_sources += " [flux from f_hot]";
        eqsys->SetOperator(id_n_re, id_f_hot, Op_nRE_fHot);
    }

    // Add source terms
    RunawaySourceTermHandler *rsth = ConstructRunawaySourceTermHandler(
        fluidGrid, hottailGrid, eqsys->GetRunawayGrid(), fluidGrid, eqsys->GetUnknownHandler(),
        eqsys->GetREFluid(), eqsys->GetIonHandler(), eqsys->GetAnalyticHottailDistribution(), oqty_terms, s
    );

    rsth->AddToOperators(Op_nRE, Op_n_tot, Op_n_i);
    desc_sources += rsth->GetDescription();

    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op_nRE, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI,
        eqsys, s, false, false,
        &oqty_terms->n_re_advective_bc, &oqty_terms->n_re_diffusive_bc,
        oqty_terms
    );
    if(hasTransport) {
        // If 'f_re' is enabled, the user should disable transport
        // for 'n_re' as it will be included automatically.
        if (transport_fre != nullptr)
            DREAM::IO::PrintWarning(
                DREAM::IO::WARNING_INCONSISTENT_RE_TRANSPORT,
                "Inconsistent runaway transport settings are used. When 'f_re' "
                "is radially transported, no transport should be applied to "
                "'n_re' as this is handled automatically."
            );
        
        desc_sources += " + transport";

        // Also enable flux limiters
        enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_r =
            (enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(MODULENAME "/adv_interp/r");
        enum OptionConstants::adv_jacobian_mode adv_jac_mode_r =
            (enum OptionConstants::adv_jacobian_mode)s->GetInteger(MODULENAME "/adv_interp/r_jac");
        real_t fluxLimiterDamping = s->GetReal(MODULENAME "/adv_interp/fluxlimiterdamping");

        Op_nRE->SetAdvectionInterpolationMethod(
            adv_interp_r, adv_jac_mode_r, FVM::FLUXGRIDTYPE_RADIAL,
            id_n_re, fluxLimiterDamping
        );
    // Account for transport of f_re...
    } else if (transport_fre != nullptr) {
        len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
        FVM::Operator *Op_fRE = new FVM::Operator(fluidGrid);
        Op_fRE->AddTerm(
            new KineticEquationTermIntegratedOverMomentum(
                fluidGrid, runawayGrid, transport_fre, id_f_re,
                eqsys->GetUnknownHandler()
            )
        );

        desc_sources += " - f_re transport";
        eqsys->SetOperator(id_n_re, id_f_re, Op_fRE);
    }

    if (!desc_sources.compare(""))
        desc_sources = "0";

    eqsys->SetOperator(id_n_re, id_n_re,  Op_nRE, "dn_re/dt =" + desc_sources);
    eqsys->SetOperator(id_n_re, id_n_tot, Op_n_tot);
    eqsys->SetOperator(id_n_re, id_n_i,   Op_n_i);

    /**
     * Load initial runaway electron density profile.
     * If the input profile is not explicitly set, then 'SetInitialValue()' is
     * called with a null-pointer which results in n_re=0 at t=0
     */
    real_t *n_re_init = LoadDataR(MODULENAME, fluidGrid->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(id_n_re, n_re_init);
    delete [] n_re_init;
}

