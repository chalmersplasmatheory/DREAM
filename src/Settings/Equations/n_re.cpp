/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Equations/TransportPrescribed.hpp"
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
    s->DefineSetting(MODULENAME "/Eceff", "Model to use for calculation of the effective critical field.", (int_t)OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL);

    DefineOptions_Transport(MODULENAME, s, false);


    s->DefineSetting(MODULENAME "/compton/mode", "Model to use for Compton seed generation.", (int_t) OptionConstants::EQTERM_COMPTON_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/compton/flux", "Gamma ray photon flux (m^-2 s^-1).", (real_t) 0.0);

    // Prescribed initial profile
    DefineDataR(MODULENAME, s, "init");

}

/**
 * Construct the equation for the runaway electron density, 'n_re'.
 * It is given by the particle flux escaping the hot-tail grid, 
 * plus any other runaway sources that are enabled.
 */
void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    len_t id_n_tot  = eqsys->GetUnknownID(OptionConstants::UQTY_N_TOT);

    // Add the transient term
    FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_nRE_2 = new FVM::Operator(fluidGrid);
    Op_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));

    std::string desc_sources = ""; 
    // Add avalanche growth rate: 
    //  - fluid mode, use analytical growth rate formula,
    //  - kinetic mode, add those knockons which are created for p>pMax 
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger(MODULENAME "/avalanche");
    // Add avalanche growth rate
    if (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID || ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID_HESSLOW){
        Op_nRE->AddTerm(new AvalancheGrowthTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );
        desc_sources += " + n_re*Gamma_ava";
    } else if ( (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) && hottailGrid ){
        // XXX: assume same momentum grid at all radii
        real_t pMax = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
        Op_nRE->AddTerm(new AvalancheSourceRP(fluidGrid, eqsys->GetUnknownHandler(),pMax, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_FLUID) );
        desc_sources += " + external avalanche";
    }
/*
AvalancheSourceRP::AvalancheSourceRP(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t pMin, RPSourceMode sm
)
*/

    // Add Dreicer runaway rate
    enum OptionConstants::eqterm_dreicer_mode dm = 
        (enum OptionConstants::eqterm_dreicer_mode)s->GetInteger(MODULENAME "/dreicer");
    switch (dm) {
        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::CONNOR_HASTIE_NOCORR, -1.0
            ));
            desc_sources += " + dreicer (CH)";
            break;

        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::CONNOR_HASTIE, -1.0
            ));
            desc_sources += " + dreicer (CH)";
            break;

        case OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK:
            Op_nRE->AddTerm(new DreicerRateTerm(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),
                eqsys->GetIonHandler(), DreicerRateTerm::NEURAL_NETWORK, -1.0
            ));
            desc_sources += " + dreicer (NN)";
            break;

        default: break;     // Don't add Dreicer runaways
    }

    // Add compton source
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)s->GetInteger(MODULENAME "/compton/mode");
    if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_FLUID){
        Op_nRE_2->AddTerm(new ComptonRateTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );
        desc_sources += " + compton";
    }
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op_nRE, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, s, false
    );
    if(hasTransport)
        desc_sources += " + transport";

    eqsys->SetOperator(id_n_re, id_n_re, Op_nRE);
    eqsys->SetOperator(id_n_re, id_n_tot, Op_nRE_2);

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
				fluidGrid, Op, id_f_hot, id_n_re, hottailGrid,
				FVM::BC::PXiExternalLoss::BOUNDARY_FLUID, bc
			));
		}

        eqsys->SetOperator(id_n_re, id_f_hot, Op_nRE_fHot, "dn_re/dt = [flux from f_hot]" + desc_sources);
    } else {
        /*FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
        Op_nRE->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqsys->SetOperator(id_n_re,id_n_re,Op_nRE, "zero");
        eqsys->initializer->AddRule(
            id_n_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );*/
    }


    /**
     * Load initial runaway electron density profile.
     * If the input profile is not explicitly set, then 'SetInitialValue()' is
     * called with a null-pointer which results in n_re=0 at t=0
     */
    real_t *n_re_init = LoadDataR(MODULENAME, fluidGrid->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(id_n_re, n_re_init);
    delete [] n_re_init;
}

