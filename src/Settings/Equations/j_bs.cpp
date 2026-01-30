/**
 * Definition of equation relating to the bootstrap current density j_bs.
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
    s->DefineSetting(MODULENAME "/init_mode", "Initialize/prescribe total current or Ohmic current", (int_t)OptionConstants::EQTERM_BOOTSTRAP_INIT_MODE_OHMIC);
}

/**
 * Construct the equation for the bootstrap current density, 'j_bs',
 * which represents j_BS / (B/Bmin) which is constant on the flux surface.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_bs(EquationSystem *eqsys, Settings *s) {

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    IonHandler *ions = eqsys->GetIonHandler();
    BootstrapCurrent *bootstrap = eqsys->GetBootstrap();
    const len_t nZ = ions->GetNZ();

    const len_t id_jbs = eqsys->GetUnknownID(OptionConstants::UQTY_J_BS);
    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_Tcold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    const len_t id_Ni = eqsys->GetUnknownID(OptionConstants::UQTY_NI_DENS);

    FVM::Operator *Op_jbs = new FVM::Operator(fluidGrid);
    Op_jbs->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.));
    eqsys->SetOperator(id_jbs, id_jbs, Op_jbs, "Redl-Sauter bootstrap");

    FVM::Operator *Op_ncold = new FVM::Operator(fluidGrid);
    Op_ncold->AddTerm( new BootstrapElectronDensityTerm(fluidGrid, unknowns, bootstrap, ions) );
    eqsys->SetOperator(id_jbs, id_ncold, Op_ncold);

    FVM::Operator *Op_Tcold = new FVM::Operator(fluidGrid);
    Op_Tcold->AddTerm( new BootstrapElectronTemperatureTerm(fluidGrid, unknowns, bootstrap, ions) );
    eqsys->SetOperator(id_jbs, id_Tcold, Op_Tcold);

    FVM::Operator *Op_Ni = new FVM::Operator(fluidGrid);
    for (len_t iZ = 0; iZ < nZ; iZ++)
        Op_Ni->AddTerm( new BootstrapIonDensityTerm(fluidGrid, unknowns, bootstrap, ions, iZ) );
    eqsys->SetOperator(id_jbs, id_Ni, Op_Ni);

    OptionConstants::uqty_T_i_eqn Ti_type = (OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if(Ti_type == OptionConstants::UQTY_T_I_INCLUDE) {
        const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);

        FVM::Operator *Op_Wi = new FVM::Operator(fluidGrid);
        for (len_t iZ = 0; iZ < nZ; iZ++)
            Op_Wi->AddTerm( new BootstrapIonThermalEnergyTerm(fluidGrid, unknowns, bootstrap, ions, iZ) );
        eqsys->SetOperator(id_jbs, id_Wi, Op_Wi);

        eqsys->initializer->AddRule(
            id_jbs,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            id_ncold, id_Tcold, id_Ni, id_Wi,
			EqsysInitializer::BOOTSTRAP
        );
    } else
        eqsys->initializer->AddRule(
            id_jbs,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            id_ncold, id_Tcold, id_Ni,
			EqsysInitializer::BOOTSTRAP
        );
}
