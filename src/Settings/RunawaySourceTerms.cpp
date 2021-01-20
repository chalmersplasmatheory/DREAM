/**
 * Initialization of runaway source terms.
 */

#include "DREAM/Equations/Fluid/TritiumRateTerm.hpp"
#include "DREAM/Equations/Fluid/ExternalAvalancheTerm.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;


RunawaySourceTermHandler *SimulationGenerator::ConstructRunawaySourceTermHandler(
    FVM::Grid *grid, FVM::Grid *hottailGrid, FVM::Grid *runawayGrid, FVM::Grid *fluidGrid,
    FVM::UnknownQuantityHandler *unknowns, RunawayFluid *REFluid,
    IonHandler *ions, Settings *s, bool signPositive
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

            if (grid == runawayGrid)
                rsth->AddSourceTerm(eqnSign + "external avalanche", new AvalancheSourceRP(grid, unknowns, pCut, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_KINETIC) );
            else if (grid == fluidGrid){
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
    if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_FLUID){
        rsth->AddSourceTerm(eqnSign + "compton", new ComptonRateTerm(grid, unknowns, REFluid, fluidGrid, -1.0) );
    }

    // Add tritium source
    bool tritium_enabled = s->GetBool(mod + "/tritium");
    if (tritium_enabled) {
        const len_t *ti = ions->GetTritiumIndices();
        for (len_t i = 0; i < ions->GetNTritiumIndices(); i++)
            rsth->AddSourceTerm(eqnSign + "tritium", new TritiumRateTerm(grid, ions, unknowns, ti[i], REFluid, -1.0));
    }

    return rsth;
}

