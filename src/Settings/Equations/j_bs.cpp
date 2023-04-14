/**
 * Definition of equations relating to the radial profile
 * of bootstrap current density j_bs.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapElectronTemperatureTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapIonDensityTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapIonThermalEnergyTerm.hpp"

#include "FVM/Equation/IdentityTerm.hpp"


using namespace DREAM;


#define MODULENAME "eqsys/j_bs"


void SimulationGenerator::DefineOptions_j_bs(Settings *s){
    s->DefineSetting(MODULENAME "/mode", "Determines which formula to use for the bootstrap current", (int_t)OptionConstants::EQTERM_BOOTSTRAP_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/bc", "Determines whether to assume plasma quantities are zero outside the edge", (int_t)OptionConstants::EQTERM_BOOTSTRAP_BC_BACKWARDS);
}

/**
 * Construct the equation for the bootstrap current density, 'j_bs',
 * which represents j_BS / (B/Bmin) which is constant on the flux surface.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_bs(
    EquationSystem *eqsys, Settings *s, struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    IonHandler *ions = eqsys->GetIonHandler();
    BootstrapCurrent *bootstrap = eqsys->GetBootstrap();

    const len_t id_jtot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
    const len_t id_jbs = eqsys->GetUnknownID(OptionConstants::UQTY_J_BS);
    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_Tcold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    const len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);


    enum OptionConstants::eqterm_bootstrap_bc bc = (OptionConstants::eqterm_bootstrap_bc)s->GetInteger("eqsys/j_bs/bc");

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.));
    eqsys->SetOperator(id_jbs, id_jbs, Op1, "Redl-Sauter bootstrap");

    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    oqty_terms->j_bs_n_cold = new BootstrapElectronDensityTerm(fluidGrid, unknowns, bootstrap, ions, bc);
    Op2->AddTerm( oqty_terms->j_bs_n_cold );
    eqsys->SetOperator(id_jbs, id_ncold, Op2);

    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
    oqty_terms->j_bs_T_cold = new BootstrapElectronTemperatureTerm(fluidGrid, unknowns, bootstrap, ions, bc);
    Op3->AddTerm( oqty_terms->j_bs_T_cold );
    eqsys->SetOperator(id_jbs, id_Tcold, Op3);

    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    oqty_terms->j_bs_N_i = new BootstrapIonDensityTerm(fluidGrid, unknowns, bootstrap, ions, bc);
    Op4->AddTerm( oqty_terms->j_bs_N_i );
    eqsys->SetOperator(id_jbs, id_Ni, Op4);

    OptionConstants::uqty_T_i_eqn Ti_type = (OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if(Ti_type == OptionConstants::UQTY_T_I_INCLUDE) {
        const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);

        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        oqty_terms->j_bs_W_i = new BootstrapIonThermalEnergyTerm(fluidGrid, unknowns, bootstrap, ions, bc);
        Op5->AddTerm( oqty_terms->j_bs_W_i );
        eqsys->SetOperator(id_jbs, id_Wi, Op5);

        eqsys->initializer->AddRule(
                id_jbs,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                id_jtot, id_ncold, id_Tcold, id_Ni, id_Wi
        );
    } else
        eqsys->initializer->AddRule(
                id_jbs,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                id_jtot, id_ncold, id_Tcold, id_Ni
        );

    // bootstrap->Rebuild();
}
