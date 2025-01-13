
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"

using namespace DREAM;
using namespace std;

#define MODULENAME "eqsys/tau_coll"


/** 
 * Implementation of an equation term representing the 
 * collision frequency nu_c = 4*pi*lnL*r0^2*c*n_cold
 */
class CollisionFrequencyTerm : public FVM::DiagonalComplexTerm {
    public:
        CoulombLogarithm *lnL;
        const len_t id_Tcold, id_ni;
        const real_t preFactor; // the constant factor nu_c / (lnL*n_cold)

        CollisionFrequencyTerm(FVM::Grid *g, FVM::UnknownQuantityHandler *u, CoulombLogarithm *lnL) 
            : FVM::DiagonalComplexTerm(g,u), lnL(lnL), 
            id_Tcold(u->GetUnknownID(OptionConstants::UQTY_T_COLD)),
            id_ni(u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES)),
            preFactor(4*M_PI*Constants::r0*Constants::r0*Constants::c)
        {
            // add unknown-quantity dependencies of the Coulomb logarithm
            AddUnknownForJacobian(u, id_Tcold);
            AddUnknownForJacobian(u, id_ni);
        }

        // set weights for this equation term
        virtual void SetWeights() override {
            for(len_t ir=0; ir<nr; ir++)
                this->weights[ir] = preFactor*lnL->evaluateAtP(ir,0);
        } 

        // derivative of coulomb logarithm with respect to T_cold and n_i
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            for(len_t n=0; n<nMultiples; n++)
                for(len_t ir=0; ir<nr; ir++)
                    this->diffWeights[n*nr + ir] = preFactor*lnL->evaluatePartialAtP(ir,0,derivId,n);
        }    
};

/**
 * Construct the equation for the time-integrated collision
 * frequency "tau_coll" appearing in the non-relativistic 
 * analytic hottail formula. This method is only called if
 * fluid hot tail is enabled with the non-relativistic distribution.
 *
 * eqsys: Equation system to put the equation in.
 */
void SimulationGenerator::ConstructEquation_tau_coll(EquationSystem *eqsys) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    const len_t id_tau = eqsys->GetUnknownID(OptionConstants::UQTY_TAU_COLL);
    const len_t id_ncold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(
        new FVM::TransientTerm(fluidGrid, id_tau, -1.0) 
    );
    Op2->AddTerm(
        new CollisionFrequencyTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid()->GetLnLambda())
    );
    eqsys->SetOperator(id_tau, id_tau, Op1);
    eqsys->SetOperator(id_tau, id_ncold, Op2, "dtau_coll/dt = nu_c");

    // Always initialize tau to 0
	const len_t nr = fluidGrid->GetNr();
	real_t *init = new real_t[nr];
	for (len_t i = 0; i < nr; i++)
		init[i] = 0.0451;

    eqsys->SetInitialValue(id_tau, init);
}
