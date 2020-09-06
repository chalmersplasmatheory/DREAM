/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/DreicerRateTerm.hpp"
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
    s->DefineSetting(MODULENAME "/avalanche", "Enable/disable secondary (avalanche) generation.", (int_t) OptionConstants::EQTERM_AVALANCHE_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/pCutAvalanche", "Minimum momentum to which the avalanche source is applied", (real_t) 0.0);
    s->DefineSetting(MODULENAME "/dreicer", "Model to use for Dreicer generation.", (int_t)OptionConstants::EQTERM_DREICER_MODE_NONE);
    s->DefineSetting(MODULENAME "/Eceff", "Model to use for calculation of the effective critical field.", (int_t)OptionConstants::COLLQTY_ECEFF_MODE_CYLINDRICAL);
    s->DefineSetting(MODULENAME "/drr", "Transport diffusion coefficient", 0, (real_t*)nullptr);

    DefineOptions_Transport(MODULENAME, s, false);

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

    // Add the transient term
    FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
    Op_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));


    // Add avalanche growth rate: 
    //  - fluid mode, use analytical growth rate formula,
    //  - kinetic mode, add those knockons which are created for p>pMax 
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger(MODULENAME "/avalanche");
    if (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID)
        Op_nRE->AddTerm(new AvalancheGrowthTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );
    else if ( (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) && hottailGrid ){
        // XXX: assume same momentum grid at all radii
        real_t pMax = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
        Op_nRE->AddTerm(new AvalancheSourceRP(fluidGrid, eqsys->GetUnknownHandler(),pMax, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_FLUID) );
    }

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

    // Prescribe transport?
    /*len_t drr_ndims[2];
    const real_t *drr = s->GetRealArray(MODULENAME "/drr", 2, drr_ndims, false);
    if (drr != nullptr) {
        s->MarkUsed(MODULENAME "/drr");

        len_t drr_nt, drr_nr;
        const real_t *drr_t = s->GetRealArray(MODULENAME "/drr_t", 1, &drr_nt);
        const real_t *drr_r = s->GetRealArray(MODULENAME "/drr_r", 1, &drr_nr);

        if (drr_ndims[0] != drr_nt || drr_ndims[1] != drr_nr)
            throw SettingsException(
                "n_re: Invalid dimensions of prescribed diffusion coefficent. Expected "
                LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT " but got "
                LEN_T_PRINTF_FMT "x" LEN_T_PRINTF_FMT ".",
                drr_nt, drr_nr, drr_ndims[0], drr_ndims[1]
            );

        // Convert 'drr' to 2D...
        const real_t **drr2d = new const real_t*[drr_nt];
        for (len_t i = 0; i < drr_nt; i++)
            drr2d[i] = drr + i*drr_nr;

        Op_nRE->AddTerm(new TransportPrescribedDiffusive(
            fluidGrid, drr_nt, drr_nr, 1, 1, drr2d, drr_t, drr_r, nullptr, nullptr,
            FVM::Interpolator3D::GRID_PXI, FVM::Interpolator3D::GRID_PXI
        ));
    }*/

    // Add transport terms, if enabled
    ConstructTransportTerm(
        Op_nRE, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, s, false
    );

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


    /**
     * Load initial runaway electron density profile.
     * If the input profile is not explicitly set, then 'SetInitialValue()' is
     * called with a null-pointer which results in n_re=0 at t=0
     */
    real_t *n_re_init = LoadDataR(MODULENAME, fluidGrid->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(id_n_re, n_re_init);
    delete [] n_re_init;
}

