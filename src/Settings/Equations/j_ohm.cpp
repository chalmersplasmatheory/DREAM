/**
 * Definition of equations relating to the radial profile 
 * of ohmic current density j_ohm.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include <algorithm>

using namespace DREAM;

/**
 * Implementation of a class which represents the sigma*E contribution to the ohmic current equation.
 */
namespace DREAM {
    class CurrentFromConductivityTerm : public FVM::WeightedIdentityTerm {
    private:
        RunawayFluid *REFluid;
        IonHandler *ionHandler;
    protected:
        virtual bool TermDependsOnUnknowns() override {return true;}
    public:
        CurrentFromConductivityTerm(FVM::Grid* g, RunawayFluid *ref, IonHandler *ih) : FVM::WeightedIdentityTerm(g), REFluid(ref), ionHandler(ih){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                //real_t w=0;
                real_t w = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir)) * REFluid->evaluateElectricalConductivity(ir, ionHandler->evaluateZeff(ir));
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}

#define MODULENAME "eqsys/j_ohm"


/**
 * Construct the equation for the ohmic current density, 'j_ohm',
 * which represents j_\Omega / (B/Bmin) which is constant on the flux surface.
 * This is zero when the hot tail grid uses collfreqmode FULL,
 * in which case the ohmic current is part of f_hot. 
 * The ohmic current j_ohm is calculated from a semi-analytical 
 * electric conductivity formula.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_j_ohm(
    EquationSystem *eqsys, Settings *s 
) {

    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();

    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");


    // If collfreqmode is FULL, the ohmic current is naturally captured in j_hot and we set j_ohm=0.
    if (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL) {
        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        //eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqn->AddTerm(new FVM::ConstantParameter(fluidGrid, 1e5));
        eqn->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));

        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn, "zero");

    // Otherwise, calculate it from conductivity
    } else {
        FVM::Equation *eqn1 = new FVM::Equation(fluidGrid);
        FVM::Equation *eqn2 = new FVM::Equation(fluidGrid);
        
        eqn2->AddTerm(new CurrentFromConductivityTerm(fluidGrid, eqsys->GetREFluid(), eqsys->GetIonHandler()));
        eqn1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_J_OHM, eqn1, "Ohmic current from conductivity");
        eqsys->SetEquation(OptionConstants::UQTY_J_OHM, OptionConstants::UQTY_E_FIELD, eqn2);
        
    }
}

