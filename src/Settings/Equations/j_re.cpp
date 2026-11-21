/**
 * Definition of equations relating to j_re (the radial profile 
 * of parallel current density j_|| / (B/B_min) of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/j_re"

/**
 * Implementation of an equation term representing the fluid runaway current
 * as a function of n_re. Here assumes that all runaways move at the  
 * speed of light in the direction of the electric field. Since n_re is the
 * flux-surface averaged density but `j_re` represents j_||re / (B/Bmin),
 * there is a geometric factor of 1/<B/Bmin> in the formula.
 */
namespace DREAM {
class RunawayFluidCurrentTerm : public FVM::DiagonalLinearTerm {
    protected:
        len_t id_Efield;
        FVM::UnknownQuantityHandler *unknowns;
        int direction;
        virtual bool TermDependsOnUnknowns() override { return direction == 0; }
        virtual void SetWeights() override{
            const real_t *Efield = direction == 0 ? unknowns->GetUnknownData(id_Efield) : nullptr;
            const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();
            for(len_t ir=0; ir<nr; ir++){
                const real_t sgn = Efield ? (Efield[ir] > 0) - (Efield[ir] < 0) : direction;
                weights[ir] = sgn * Constants::ec * Constants::c / FSA_B[ir];
            }
        }
    public:
    RunawayFluidCurrentTerm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, int direction=0)
        : FVM::DiagonalLinearTerm(g), unknowns(u), direction(direction) {
        id_Efield = u->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    }
};
}

/**
 * Construct the equation for the runaway parallel current, 'j_re'.
 * Include the moment of the resolved runaway distribution, when present.
 * Fluid runaways contribute a current at c, using their stored pitch sign
 * with the two-way boundary and the electric-field sign otherwise.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_re(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    len_t id_j_re = eqsys->GetUnknownID(OptionConstants::UQTY_J_RE);
    len_t id_n_re = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);


    // Identity part
    FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
    eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    eqsys->SetOperator(id_j_re, id_j_re, eqnIdent);

    auto *unknowns = eqsys->GetUnknownHandler();
    len_t id_f_re = 0;

    if (runawayGrid) {
        id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
        auto *moment = new FVM::Operator(fluidGrid);
        moment->AddTerm(new CurrentDensityFromDistributionFunction(
            fluidGrid, runawayGrid, id_j_re, id_f_re, unknowns
        ));
        eqsys->SetOperator(id_j_re, id_f_re, moment, "Moment of f_re");
    }

    if (unknowns->HasUnknown(OptionConstants::UQTY_N_EXTERNAL_POS)) {
        const len_t idPositive = eqsys->GetUnknownID(OptionConstants::UQTY_N_EXTERNAL_POS);
        const len_t idNegative = eqsys->GetUnknownID(OptionConstants::UQTY_N_EXTERNAL_NEG);
        auto *positive = new FVM::Operator(fluidGrid);
        positive->AddTerm(new RunawayFluidCurrentTerm(fluidGrid, unknowns, +1));
        eqsys->SetOperator(id_j_re, idPositive, positive, runawayGrid
            ? "Moment of f_re plus signed external current" : "Signed external current");
        auto *negative = new FVM::Operator(fluidGrid);
        negative->AddTerm(new RunawayFluidCurrentTerm(fluidGrid, unknowns, -1));
        eqsys->SetOperator(id_j_re, idNegative, negative);
        if (runawayGrid)
            eqsys->initializer->AddRule(id_j_re, EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr, id_f_re, idPositive, idNegative);
        else
            eqsys->initializer->AddRule(id_j_re, EqsysInitializer::INITRULE_EVAL_EQUATION,
                nullptr, idPositive, idNegative);
    } else if (runawayGrid) {
        eqsys->initializer->AddRule(id_j_re, EqsysInitializer::INITRULE_EVAL_EQUATION, nullptr, id_f_re);
    } else {
        auto *eqn = new FVM::Operator(fluidGrid);
        eqn->AddTerm(new RunawayFluidCurrentTerm(fluidGrid, unknowns));
        eqsys->SetOperator(id_j_re, id_n_re, eqn, "j_re = sgn(E)*e*c*n_re");
        eqsys->initializer->AddRule(id_j_re, EqsysInitializer::INITRULE_EVAL_EQUATION, nullptr, id_n_re);
    }
}
