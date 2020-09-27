/**
 * Definition of equations relating to n_cold.
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/FreeElectronDensityTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/n_cold"

/**
 * Construct the equation for the cold electron density, 'n_cold'.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_n_cold(
    EquationSystem *eqsys, Settings *s
) {
    enum OptionConstants::uqty_n_cold_eqn eqn = (enum OptionConstants::uqty_n_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (eqn) {
        case OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED:
            ConstructEquation_n_cold_prescribed(eqsys, s);
            break;
        
        case OptionConstants::UQTY_N_COLD_EQN_SELFCONSISTENT:
            ConstructEquation_n_cold_selfconsistent(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for 'n_cold': %d.",
                eqn
            );
    }
}

/**
 * Construct the equation describing a prescribed cold
 * electron density.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct
 *        the equation.
 */
void SimulationGenerator::ConstructEquation_n_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Operator *eqn = new FVM::Operator(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT_intp(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetOperator(OptionConstants::UQTY_N_COLD, OptionConstants::UQTY_N_COLD, eqn, "Prescribed");
    //eqsys->SetInitialValue(OptionConstants::UQTY_N_COLD, interp->Eval(t0), t0);

    eqsys->initializer->AddRule(
        OptionConstants::UQTY_N_COLD,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}

/**
 * Construct the equation describing the cold electron density as
 * (if superthermal grid and hot electrons are transferred to n_cold via lower p boundary:)
 *   n_cold = n_free - n_hot - n_re
 * (if particle conserving collision operator:)
 *   n_cold = n_hot + n_re
 */
void SimulationGenerator::ConstructEquation_n_cold_selfconsistent(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_nhot  = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    const len_t id_nre   = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    const len_t id_ni    = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);


    if (eqsys->GetHotTailGrid()){
        enum OptionConstants::collqty_collfreq_mode collfreq_mode =
            (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");

        // If conservative hot tail grid, n_cold and n_hot are the same quantity 
        // (since n_hot is the thermal population in that case) 
        if (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
            FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
            FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
            FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
            Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
            Op2->AddTerm(new FVM::IdentityTerm(fluidGrid));
            Op3->AddTerm(new FVM::IdentityTerm(fluidGrid));

            eqsys->SetOperator(id_ncold, id_ncold, Op1, "n_cold = n_hot + n_re");
            eqsys->SetOperator(id_ncold, id_nhot, Op2);
            eqsys->SetOperator(id_ncold, id_nre, Op3);
            
            // Initialization
            eqsys->initializer->AddRule(
                OptionConstants::UQTY_N_COLD,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_nhot,
                id_nre
            );

        
        // Otherwise, n_cold are the thermal particles who are not fast or RE 
        } else {
            FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
            FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
            FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
            FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
            
            Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
            Op2->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
            Op3->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
            Op4->AddTerm(new FreeElectronDensityTerm(fluidGrid,eqsys->GetIonHandler()));
            eqsys->SetOperator(id_ncold, id_ncold, Op1, "n_cold = n_free - n_hot - n_re");
            eqsys->SetOperator(id_ncold, id_nhot, Op2);
            eqsys->SetOperator(id_ncold, id_nre, Op3);
            eqsys->SetOperator(id_ncold, id_ni, Op4);

            // Initialization
            eqsys->initializer->AddRule(
                OptionConstants::UQTY_N_COLD,
                EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr,
                // Dependencies
                id_ni, id_nre, id_nhot
            );
        }
    
    // If no hot tail grid, n_cold is the free electron density minus runaways
    } else {
        FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
        Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        Op2->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        Op3->AddTerm(new FreeElectronDensityTerm(fluidGrid,eqsys->GetIonHandler()));
        eqsys->SetOperator(id_ncold, id_ncold, Op1, "n_cold = n_free - n_re");
        eqsys->SetOperator(id_ncold, id_nre, Op2);
        eqsys->SetOperator(id_ncold, id_ni, Op3);

        // Initialization
        eqsys->initializer->AddRule(
            OptionConstants::UQTY_N_COLD,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            // Dependencies
            id_ni, id_nre
        );
    }


}

