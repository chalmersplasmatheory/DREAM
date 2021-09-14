/**
 * Definition of equations relating to W_hot (the kinetic energy density of hot electrons).
 */

#include <sstream>
#include <iomanip>

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/Fluid/KineticEnergyFromDistributionFunction.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


/**
 * Construct the equation for the hot kinetic energy density, 'W_hot'.
 *
 * eqsys:  Equation system to put the equation in.
 * s:      Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_W_hot(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid   = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    len_t id_W_hot = eqsys->GetUnknownID(OptionConstants::UQTY_W_HOT);
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    // Identity part
    FVM::Operator *Op0 = new FVM::Operator(fluidGrid);
    Op0->AddTerm(new FVM::IdentityTerm(fluidGrid, -1.0));
    eqsys->SetOperator(id_W_hot, id_W_hot, Op0);

    std::string desc = "Energy moment of f_hot";

    FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
    real_t pThreshold = 0.0;
    enum OptionConstants::collqty_collfreq_mode collfreq_mode =
        (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
        // With collfreq_mode FULL, n_hot is defined as density above some threshold.
        pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        
        std::ostringstream str;
        str <<std::fixed << std::setprecision(3) << pThreshold;
        switch(pMode){
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p>"+str.str()+"*me*c)";
                break;
            } 
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p>"+str.str()+"*pThermal)";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL_SMOOTH:{
                desc = "integral(me*c^2(gamma-1)*f_hot), smooth lower limit at p="+str.str()+"*pThermal";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_MC:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p<"+str.str()+"*me*c)";
                break;
            } 
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL:{
                desc = "integral(me*c^2(gamma-1)*f_hot, p<"+str.str()+"*pThermal)";
                break;
            }
            case FVM::MomentQuantity::P_THRESHOLD_MODE_MAX_THERMAL_SMOOTH:{
                desc = "integral(me*c^2(gamma-1)*f_hot), smooth upper limit at p="+str.str()+"*pThermal";
                break;
            }    
        }
    }
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new KineticEnergyFromDistributionFunction(
            fluidGrid, hottailGrid, id_W_hot, id_f_hot,eqsys->GetUnknownHandler(),pThreshold, pMode)
        );
    eqsys->SetOperator(id_W_hot, id_f_hot, Op1, desc);

    // Set initialization method
    eqsys->initializer->AddRule(
        id_W_hot,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        // Dependencies
        id_f_hot
    );
}

