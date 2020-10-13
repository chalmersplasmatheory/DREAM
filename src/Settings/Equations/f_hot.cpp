/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

#include <iostream>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Kinetic/BCIsotropicSourcePXi.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"
#include "DREAM/Equations/Kinetic/PitchScatterTerm.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "DREAM/Equations/Kinetic/ParticleSourceTerm.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;
using namespace std;

#define MODULENAME "eqsys/f_hot"


/**
 * Define settings for the hot-tail distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_hot(Settings *s) {
    DefineOptions_f_general(s, MODULENAME);

    // Cold electron definition
    s->DefineSetting(MODULENAME "/pThreshold", "Threshold momentum that defines n_hot from f_hot when resolving thermal population on grid.", (real_t) 10.0);
    s->DefineSetting(MODULENAME "/pThresholdMode", "Unit of provided threshold momentum pThreshold (thermal or mc).", (int_t) FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL);

}

/**
 * Construct the equation for the hot electron distribution
 * function. This function is only called if the hot-tail grid
 * is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_f_hot(
    EquationSystem *eqsys, Settings *s
) {
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
    bool addExternalBC = not eqsys->HasRunawayGrid();
    bool addInternalBC = true;
    FVM::Operator *eqn = ConstructEquation_f_general(
        s, MODULENAME, eqsys, id_f_hot, hottailGrid, eqsys->GetHotTailGridType(),
        eqsys->GetHotTailCollisionHandler(), addExternalBC, addInternalBC
    );

    
    // Add kinetic-kinetic boundary condition if necessary...
    if (!addExternalBC) {
		len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			hottailGrid, hottailGrid, eqsys->GetRunawayGrid(), eqn,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
		));
    }

    // PARTICLE SOURCE TERMS
    const len_t id_Sp = eqsys->GetUnknownID(OptionConstants::UQTY_S_PARTICLE);
    const len_t id_n_cold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_n_hot  = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    FVM::Operator *Op_source = new FVM::Operator(hottailGrid);
    ParticleSourceTerm::ParticleSourceShape sourceShape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_MAXWELLIAN;
//    ParticleSourceTerm::ParticleSourceShape sourceShape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_DELTA;
    Op_source->AddTerm(new ParticleSourceTerm(hottailGrid,eqsys->GetUnknownHandler(),sourceShape) );
    eqsys->SetOperator(id_f_hot, id_Sp, Op_source);

    // Temporarily switch the self-consistent particle source by prescribing it to 0:
    bool useParticleSourceTerm = true; 
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    OptionConstants::collqty_collfreq_mode collfreq_mode = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(useParticleSourceTerm && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL)){
        // Set self-consistent particle source equation: Require total f_hot density to equal n_cold+n_hot
        FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
        FVM::Operator *Op4 = new FVM::Operator(fluidGrid);

        Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        Op2->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
        Op3->AddTerm(new DensityFromDistributionFunction(
                fluidGrid, hottailGrid, id_Sp, id_f_hot,eqsys->GetUnknownHandler())
            );
        // Add regularizing term that prevents system from becoming singular when density is naturally conserved
        const real_t REGULARIZATION_FACTOR = 1e-5;
        Op4->AddTerm(new FVM::IdentityTerm(fluidGrid,-REGULARIZATION_FACTOR));
        eqsys->SetOperator(id_Sp, id_n_cold, Op1, "integral(f_hot) = n_cold + n_hot + adhoc*S_p");
        eqsys->SetOperator(id_Sp, id_n_hot,  Op2);
        eqsys->SetOperator(id_Sp, id_f_hot,  Op3);
        eqsys->SetOperator(id_Sp, id_Sp,  Op4);

    // if inactivated, just set to 0
    } else {
        FVM::Operator *Op = new FVM::Operator(fluidGrid);
        Op->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqsys->SetOperator(id_Sp, id_Sp, Op, "zero");
    }
    // initialize Sp to 0
    real_t *initZero = new real_t[fluidGrid->GetNr()];
    for(len_t ir=0; ir<fluidGrid->GetNr();ir++)
        initZero[ir] = 0;
    eqsys->SetInitialValue(id_Sp, initZero);
    delete [] initZero;   
}

