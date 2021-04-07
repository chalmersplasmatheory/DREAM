/**
 * Definition of equations relating to j_re (the radial profile 
 * of parallel current density j_|| / (B/B_min) of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/CurrentDensityFromAnalyticRE.hpp"
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
        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void SetWeights() override{
            const real_t *Efield = unknowns->GetUnknownData(id_Efield);
            const real_t *FSA_B = this->grid->GetRadialGrid()->GetFSA_B();
            for(len_t ir=0; ir<nr; ir++){
                real_t sgn = (Efield[ir] > 0) - (Efield[ir] < 0);
                weights[ir] = sgn * Constants::ec * Constants::c / FSA_B[ir];
            }
        }
    public:
     RunawayFluidCurrentTerm(FVM::Grid* g,FVM::UnknownQuantityHandler *u)
    : FVM::DiagonalLinearTerm(g), unknowns(u) {
        id_Efield = u->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    }
};
}

/**
 * Construct the equation for the runaway parallel current, 'j_re'.
 * If the runaway grid is enabled, j_re will be an integral of
 * the runaway electron distribution. If it does not exist, it is set
 * to e*c*n_re.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_re(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    len_t id_j_re = eqsys->GetUnknownID(OptionConstants::UQTY_J_RE);
    len_t id_n_re = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);


    // Identity part
    FVM::Operator *eqnIdent = new FVM::Operator(fluidGrid);
    eqnIdent->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    eqsys->SetOperator(id_j_re, id_j_re, eqnIdent);

    FVM::Operator *Op = new FVM::Operator(fluidGrid);

    OptionConstants::uqty_distribution_mode REDistMode = (OptionConstants::uqty_distribution_mode) s->GetInteger("eqsys/f_re/mode");
    // if runawayGrid is enabled, take moment of f_re, otherwise e*c*n_re
    // TODO: if integral(f_re) significantly deviates from n_re, warn that
    // the RE current is not well resolved?
    if (runawayGrid) {
        len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
        Op->AddTerm(new CurrentDensityFromDistributionFunction(
            fluidGrid, runawayGrid, id_j_re, id_f_re, eqsys->GetUnknownHandler()
        ));
        eqsys->SetOperator(id_j_re, id_f_re, Op, "e*v_|| moment of f_re");

        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            // Dependencies
            id_f_re
        );

    } else {
        std::string desc;
        if (REDistMode == OptionConstants::UQTY_DISTRIBUTION_MODE_ANALYTICAL){
            Op->AddTerm(new CurrentDensityFromAnalyticRE(
                fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetAnalyticREDistribution()
            ));
            desc = "j_re = e*v_|| moment of analytic f_re";
        } else { 
            Op->AddTerm(new RunawayFluidCurrentTerm(
                fluidGrid, eqsys->GetUnknownHandler()
            ));
            desc = "j_re = sgn(E)*e*c*n_re";
        }
        eqsys->SetOperator(id_j_re, id_n_re, Op, desc);
        // Set initialization method
        eqsys->initializer->AddRule(
            id_j_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION,
            nullptr,
            id_n_re
        );
    }

}
