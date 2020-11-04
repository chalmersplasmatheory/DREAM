/**
 * Initialization of runaway source terms.
 */

#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;


RunawaySourceTermHandler *SimulationGenerator::ConstructRunawaySourceTermHandler(
    FVM::Grid *grid, FVM::Grid *hottailGrid, FVM::UnknownQuantityHandler *unknowns,
    RunawayFluid *REFluid, IonHandler *ions, Settings *s
) {
    const std::string &mod = "eqsys/n_re";

    RunawaySourceTermHandler *rsth = new RunawaySourceTermHandler();

    // Add avalanche growth rate: 
    //  - fluid mode, use analytical growth rate formula,
    //  - kinetic mode, add those knockons which are created for p>pMax 
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger(mod + "/avalanche");
    // Add avalanche growth rate
    if (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID || ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID_HESSLOW){
        rsth->AddSourceTerm(" + n_re*Gamma_ava", new AvalancheGrowthTerm(grid, unknowns, REFluid,-1.0) );
    } else if ( (ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) && hottailGrid ){
        // XXX: assume same momentum grid at all radii
        real_t pMax = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
        rsth->AddSourceTerm(" + external avalanche", new AvalancheSourceRP(grid, unknowns, pMax, -1.0, AvalancheSourceRP::RP_SOURCE_MODE_FLUID) );
    }

    // Add Dreicer runaway rate
    enum OptionConstants::eqterm_dreicer_mode dm = 
        (enum OptionConstants::eqterm_dreicer_mode)s->GetInteger(mod + "/dreicer");
    switch (dm) {
        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR:
            rsth->AddSourceTerm(" + dreicer (CH)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::CONNOR_HASTIE_NOCORR, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_CONNOR_HASTIE:
            rsth->AddSourceTerm(" + dreicer (CH)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::CONNOR_HASTIE, -1.0
            ));
            break;

        case OptionConstants::EQTERM_DREICER_MODE_NEURAL_NETWORK:
            rsth->AddSourceTerm(" + dreicer (NN)", new DreicerRateTerm(
                grid, unknowns, REFluid,
                ions, DreicerRateTerm::NEURAL_NETWORK, -1.0
            ));
            break;

        default: break;     // Don't add Dreicer runaways
    }

    // Add compton source
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)s->GetInteger(mod + "/compton/mode");
    if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_FLUID){
        rsth->AddSourceTerm(" + compton", new ComptonRateTerm(grid, unknowns, REFluid, -1.0) );
    }

    return rsth;
}

