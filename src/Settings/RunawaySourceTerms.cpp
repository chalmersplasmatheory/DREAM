/**
 * Initialization of runaway source terms.
 */

#include "DREAM/Equations/Fluid/TritiumRateTerm.hpp"
#include "DREAM/Equations/Fluid/HottailRateTermHighZ.hpp"
#include "DREAM/Equations/Fluid/LCFSLossRateTerm.hpp"
#include "DREAM/Equations/Fluid/ExternalAvalancheTerm.hpp"
#include "DREAM/Equations/Kinetic/ComptonSource.hpp"
#include "DREAM/Equations/Kinetic/TritiumSource.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;


RunawaySourceTermHandler *SimulationGenerator::ConstructRunawaySourceTermHandler(
    FVM::Grid *grid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid, FVM::Grid *fluidGrid,
    FVM::UnknownQuantityHandler *unknowns, RunawayFluid *REFluid,
    IonHandler *ions, AnalyticDistributionHottail *distHT, 
    struct OtherQuantityHandler::eqn_terms *oqty_terms, Settings *s, bool signPositive
) {
    const std::string &mod = "eqsys/n_re";

    std::string eqnSign = signPositive ? " + " : " - ";
    
    

    RunawaySourceTermHandler *rsth = new RunawaySourceTermHandler();

    // Add avalanche growth rate: 
    //  - fluid mode, use analytical growth rate formula,
    //  - kinetic mode, add those knockons which are created for p>pMax 
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger(mod + "/avalanche");
    // Add avalanche growth rate
    if (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID || ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID_HESSLOW)
        rsth->AddSourceTerm(eqnSign + "n_re*Gamma_ava", new AvalancheGrowthTerm(grid, unknowns, REFluid, fluidGrid, -1.0) );

    else if (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) {
        if (hottailGrid || runawayGrid != nullptr) {
            // XXX: assume same momentum grid at all radii
            real_t pCut;
            if (hottailGrid != nullptr)
                pCut = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
            else if (runawayGrid != nullptr)
                pCut = runawayGrid->GetMomentumGrid(0)->GetP1_f(0);

                if (grid == runawayGrid) {
                    rsth->AddSourceTerm(eqnSign + "external avalanche", new AvalancheSourceRP(grid, unknowns, pCut, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_KINETIC, AvalancheSourceRP::RP_SOURCE_PITCH_POSITIVE) );
                    rsth->AddAvalancheNreNeg(new AvalancheSourceRP(grid, unknowns, pCut, +1.0, AvalancheSourceRP::RP_SOURCE_MODE_KINETIC, AvalancheSourceRP::RP_SOURCE_PITCH_POSITIVE));
                    rsth->AddAvalancheNreNegPos(new AvalancheSourceRP(grid, unknowns, pCut, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_KINETIC, AvalancheSourceRP::RP_SOURCE_PITCH_NEGATIVE));
                } else if (grid == fluidGrid){
                    if(runawayGrid == nullptr) // match external growth to fluid formula in E~<Eceff limit
                        rsth->AddSourceTerm(eqnSign + "external avalanche", new ExternalAvalancheTerm(grid, pCut, -2.0, REFluid, unknowns, -1.0)  );
                    else  // use regular external RE growth (RP integrated over p>pCut)
                        rsth->AddSourceTerm(eqnSign + "external avalanche", new AvalancheSourceRP(grid, unknowns, pCut, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_FLUID) );
                }
        } else
            DREAM::IO::PrintWarning(DREAM::IO::WARNING_KINETIC_AVALANCHE_NO_HOT_GRID, "A kinetic avalanche term is used, but the hot-tail grid is disabled. Ignoring avalanche source...");
    }

    // Add Dreicer runaway rate
    enum OptionConstants::eqterm_dreicer_mode dm = 
        (enum OptionConstants::eqterm_dreicer_mode)s->GetInteger(mod + "/dreicer");
    switch (dm) {
        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR:
            rsth->AddSourceTerm(eqnSign + "dreicer (CH)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::CONNOR_HASTIE_NOCORR, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE:
            rsth->AddSourceTerm(eqnSign + "dreicer (CH)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::CONNOR_HASTIE, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK:
            rsth->AddSourceTerm(eqnSign + "dreicer (NN)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::NEURAL_NETWORK, -1.0
            ));
            break;

        default: break;     // Don't add Dreicer runaways
    }

    // Add compton source
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)s->GetInteger(mod + "/compton/mode");
    if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_FLUID)
        rsth->AddSourceTerm(eqnSign + "compton", new ComptonRateTerm(grid, unknowns, REFluid, fluidGrid, -1.0) );
    else if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_KINETIC) {
        if (hottailGrid || runawayGrid != nullptr) {
            real_t pLower;
            if (hottailGrid != nullptr)
                pLower = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
                
            else if (runawayGrid != nullptr)
                pLower = runawayGrid->GetMomentumGrid(0)->GetP1_f(0);
            
            if(grid == fluidGrid) {
                oqty_terms->comptonSource_fluid = new ComptonSource(grid, unknowns, LoadDataT("eqsys/n_re/compton", s, "flux"),
                    s->GetReal("eqsys/n_re/compton/gammaInt"), s->GetReal("eqsys/n_re/compton/C1"), s->GetReal("eqsys/n_re/compton/C2"), s->GetReal("eqsys/n_re/compton/C3"), 
                    pLower, -1.0, ComptonSource::SOURCE_MODE_FLUID, REFluid);
                rsth->AddSourceTerm(eqnSign + "fluid Compton", oqty_terms->comptonSource_fluid);
            } else {
                ComptonSource *comptonSource = new ComptonSource(grid, unknowns, LoadDataT("eqsys/n_re/compton", s, "flux"), 
                    s->GetReal("eqsys/n_re/compton/gammaInt"), s->GetReal("eqsys/n_re/compton/C1"), s->GetReal("eqsys/n_re/compton/C2"), s->GetReal("eqsys/n_re/compton/C3"), 
                    pLower, -1.0, ComptonSource::SOURCE_MODE_KINETIC);
                rsth->AddSourceTerm(eqnSign + "kinetic Compton", comptonSource);
                if (grid == runawayGrid) 
                    oqty_terms->comptonSource_runaway = comptonSource;
            }
        } else {
            DREAM::IO::PrintWarning(DREAM::IO::WARNING_KINETIC_COMPTON_NO_HOT_GRID, "A kinetic Compton term is used, but the hot-tail grid and runaway grids are disabled. Ignoring Compton source...");
        }
    }
    
    // Add tritium source: 
    OptionConstants::eqterm_tritium_mode tritium_mode = (enum OptionConstants::eqterm_tritium_mode)s->GetInteger(mod + "/tritium");
    if (tritium_mode == OptionConstants::EQTERM_TRITIUM_MODE_FLUID){
        const len_t *ti = ions->GetTritiumIndices();
        for (len_t i = 0; i < ions->GetNTritiumIndices(); i++)
            rsth->AddSourceTerm(eqnSign + "fluid tritium", new TritiumRateTerm(grid, ions, unknowns, ti[i], REFluid, -1.0));
    } else if (tritium_mode == OptionConstants::EQTERM_TRITIUM_MODE_KINETIC) {
        if (hottailGrid || runawayGrid != nullptr) {
            const len_t *ti = ions->GetTritiumIndices();
            real_t pLower;
            if (hottailGrid != nullptr)
                pLower = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
            else if (runawayGrid != nullptr)
                pLower = runawayGrid->GetMomentumGrid(0)->GetP1_f(0);
            
            if(grid == fluidGrid) {
                for (len_t iT = 0; iT < ions->GetNTritiumIndices(); iT++){
                    rsth->AddSourceTerm(eqnSign + "kinetic tritium", new TritiumSource(grid, unknowns, ions, ti[iT], pLower, -1.0, TritiumSource::SOURCE_MODE_FLUID));
                }
            } else {
                for (len_t iT = 0; iT < ions->GetNTritiumIndices(); iT++){
                    rsth->AddSourceTerm(eqnSign + "kinetic tritium", new TritiumSource(grid, unknowns, ions, ti[iT], pLower, -1.0, TritiumSource::SOURCE_MODE_KINETIC));
                }
            }
        } else {
            DREAM::IO::PrintWarning(DREAM::IO::WARNING_KINETIC_TRITIUM_NO_HOT_GRID, "A kinetic tritium term is used, but the hot-tail and runaway grids are disabled. Ignoring tritium source...");
        }
    } 

    // Add hottail source
    OptionConstants::eqterm_hottail_mode hottail_mode = (enum OptionConstants::eqterm_hottail_mode)s->GetInteger(mod + "/hottail");
    if(distHT!=nullptr && hottail_mode == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC){
        oqty_terms->n_re_hottail_rate = new HottailRateTermHighZ(
            grid, distHT, unknowns, ions, 
            REFluid->GetLnLambda(), -1.0
        );
        rsth->AddSourceTerm(eqnSign + "hottail", oqty_terms->n_re_hottail_rate);
    } else if (distHT!=nullptr && hottail_mode == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC){
        throw SettingsException(
			"The 'HOTTAIL_MODE_ANALYTIC' has not been implemented yet."
		);
    }
     
    // Add LCFS loss term (add exception if (nr==1 and lcfs_mode != LCFS_LOSS_MODE_DISABLED)...)
    OptionConstants::eqterm_lcfs_loss_mode lcfs_mode = (enum OptionConstants::eqterm_lcfs_loss_mode)s->GetInteger(mod + "/lcfs_loss");
	if (lcfs_mode != OptionConstants::EQTERM_LCFS_LOSS_MODE_DISABLED) {
		if (grid->GetNr() == 1){
			throw SettingsException(
				"The LCFS loss term is not compatible with using only one radial grid point."
			);
		}

		bool lcfs_user_input_psi = (len_t)s->GetInteger(mod + "/lcfs_user_input_psi");
		real_t lcfs_psi_edge_t0 = s->GetReal(mod + "/lcfs_psi_edge_t0");
		real_t *lcfs_t_loss = LoadDataR("eqsys/n_re", grid->GetRadialGrid(), s, "lcfs_t_loss");

		if(lcfs_mode == OptionConstants::EQTERM_LCFS_LOSS_MODE_FLUID){
			oqty_terms->lcfsLossRate_fluid = new LCFSLossRateTerm(
				grid, unknowns, fluidGrid, -1.0, lcfs_t_loss,
				lcfs_user_input_psi, lcfs_psi_edge_t0
			);
			rsth->AddSourceTerm(eqnSign + "n_re*lcfs_loss", oqty_terms->lcfsLossRate_fluid);
		} else if (lcfs_mode == OptionConstants::EQTERM_LCFS_LOSS_MODE_KINETIC){
			rsth->AddSourceTerm(
				eqnSign + "n_re*lcfs_loss",
				new LCFSLossRateTerm(
					grid, unknowns, runawayGrid, -1.0, lcfs_t_loss,
					lcfs_user_input_psi, lcfs_psi_edge_t0
				)
			);
		}
	}
    
    
    return rsth;
}

