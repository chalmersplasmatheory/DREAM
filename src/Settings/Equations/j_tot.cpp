/**
 * Definition of equations relating to j_tot, 
 * representing the total plasma current:
 *  j_tot = j_ohm + j_hot + j_RE
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_tot"

/**
 * Define options for the electric field module.
 */
void SimulationGenerator::DefineOptions_j_tot(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED);
}



/**
 * Construct the equation for the total plasma current density, 'j_tot',
 * and for the total plasma current 'I_p'.
 *
 * TODO: Proper equation for j_RE (integral of f_RE + e*c*n_RE_external)
 * 
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_j_tot(
    EquationSystem *eqsys, Settings* /* s */
) {
    const len_t id_j_tot = eqsys->GetUnknownID(OptionConstants::UQTY_J_TOT);
    const len_t id_j_ohm = eqsys->GetUnknownID(OptionConstants::UQTY_J_OHM);
    const len_t id_j_hot = eqsys->GetUnknownID(OptionConstants::UQTY_J_HOT);
    const len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    const len_t id_I_p   = eqsys->GetUnknownID(OptionConstants::UQTY_I_P);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    FVM::Operator *eqn0 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn1 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn2 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn3 = new FVM::Operator(fluidGrid);

    eqn0->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn2->AddTerm(new FVM::IdentityTerm(fluidGrid));
    eqn3->AddTerm(new FVM::IdentityTerm(fluidGrid, Constants::ec * Constants::c));
    
    eqsys->SetOperator(id_j_tot, id_j_tot, eqn0, "j_tot = j_ohm + j_hot + e*c*n_RE");
    eqsys->SetOperator(id_j_tot, id_j_ohm, eqn1);
    eqsys->SetOperator(id_j_tot, id_j_hot, eqn2);
    eqsys->SetOperator(id_j_tot, id_n_re, eqn3);

    // Initialization
    eqsys->initializer->AddRule(
        id_j_tot,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_j_ohm, id_j_hot, id_n_re
    );



    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    /**
     * Set equation for the total plasma current I_p
     * (as an integral over j_tot).
     */
    FVM::Operator *eqn_Ip1 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Ip2 = new FVM::Operator(scalarGrid);
    
    eqn_Ip1->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqn_Ip2->AddTerm(new TotalPlasmaCurrentFromJTot(scalarGrid,fluidGrid,eqsys->GetUnknownHandler(),id_j_tot));

    eqsys->SetOperator(id_I_p, id_I_p,   eqn_Ip1, "Ip = integral(j_tot)");
    eqsys->SetOperator(id_I_p, id_j_tot, eqn_Ip2);

    eqsys->initializer->AddRule(
        id_I_p,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_j_tot
    );

}


