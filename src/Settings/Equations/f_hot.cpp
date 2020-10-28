/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

#include <iostream>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "DREAM/Equations/Fluid/FreeElectronDensityTransientTerm.hpp"
#include "DREAM/Equations/Fluid/KineticEquationTermIntegratedOverMomentum.hpp"
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

// different methods for modelling the particle source S_particle
enum particleSourceType {
    PARTICLE_SOURCE_ZERO     = 1,
    PARTICLE_SOURCE_IMPLICIT = 2,
    PARTICLE_SOURCE_EXPLICIT = 3
};

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
    s->DefineSetting(MODULENAME "/particleSource", "Include particle source which enforces the integral over the distribution to follow n_hot+n_cold.", (int_t) PARTICLE_SOURCE_EXPLICIT);

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
    EquationSystem *eqsys, Settings *s, struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
    bool addExternalBC = not eqsys->HasRunawayGrid();
    bool addInternalBC = true;
    bool rescaleMaxwellian = true;
    FVM::Operator *eqn = ConstructEquation_f_general(
        s, MODULENAME, eqsys, id_f_hot, hottailGrid, eqsys->GetHotTailGridType(),
        eqsys->GetHotTailCollisionHandler(), addExternalBC, addInternalBC,
        &oqty_terms->f_hot_advective_bc, &oqty_terms->f_hot_diffusive_bc,rescaleMaxwellian
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
    FVM::Operator *Op_source = new FVM::Operator(hottailGrid);
    ParticleSourceTerm::ParticleSourceShape sourceShape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_MAXWELLIAN;
//    ParticleSourceTerm::ParticleSourceShape sourceShape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_DELTA;
    Op_source->AddTerm(new ParticleSourceTerm(hottailGrid,eqsys->GetUnknownHandler(),sourceShape) );
    eqsys->SetOperator(id_f_hot, id_Sp, Op_source);

    // Enable particle source term ?
    particleSourceType particleSource = (particleSourceType) s->GetInteger(MODULENAME "/particleSource"); 
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    OptionConstants::collqty_collfreq_mode collfreq_mode = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(particleSource==PARTICLE_SOURCE_IMPLICIT && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL))
        ConstructEquation_S_particle_implicit(eqsys, s);
    else if(particleSource==PARTICLE_SOURCE_EXPLICIT && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL))
        ConstructEquation_S_particle_explicit(eqsys, s, oqty_terms);
    else {
        // if inactivated, just prescribe to 0
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

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Rosenbluth-Putvinski source
 */
namespace DREAM {
    class TotalElectronDensityFromKineticAvalanche : public FVM::DiagonalQuadraticTerm {
    public:
        real_t pLower, scaleFactor;
        TotalElectronDensityFromKineticAvalanche(FVM::Grid* g, real_t pLower, FVM::UnknownQuantityHandler *u, real_t scaleFactor = 1.0) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_TOT),u), pLower(pLower), scaleFactor(scaleFactor) {}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = scaleFactor * AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(grid->GetRadialGrid()->GetFSA_B(i),pLower);
        }
    };
}

/**
 * Build the equation for S_particle, which contains the rate at which the total local free electron density changes
 * (i.e. by ionization, transport of hot electrons, runaway sources).
 * Uses the implicit method, where the source amplitude is indirectly set by requiring
 * that the hot distribution integrates to n_cold + n_hot.
 */
void SimulationGenerator::ConstructEquation_S_particle_implicit(EquationSystem *eqsys, Settings *){
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    const len_t id_Sp = eqsys->GetUnknownID(OptionConstants::UQTY_S_PARTICLE);
    const len_t id_n_cold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const len_t id_n_hot = eqsys->GetUnknownID(OptionConstants::UQTY_N_HOT);
    const len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    Op2->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));
    Op3->AddTerm(new DensityFromDistributionFunction(
            fluidGrid, eqsys->GetHotTailGrid(), id_Sp, id_f_hot,eqsys->GetUnknownHandler()
        ));
    eqsys->SetOperator(id_Sp, id_n_cold, Op1, "integral(f_hot) = n_cold + n_hot");
    eqsys->SetOperator(id_Sp, id_n_hot,  Op2);
    eqsys->SetOperator(id_Sp, id_f_hot,  Op3);
}

/**
 * Build the equation for S_particle, which contains the rate at which the total local free electron density changes
 * (i.e. by ionization, transport of hot electrons, runaway sources).
 * Uses the explicit model, where the required amplitude of S_particle is directly set
 * by setting it to the sum of all equation terms that modify the electron density
 */
void SimulationGenerator::ConstructEquation_S_particle_explicit(EquationSystem *eqsys, Settings *s, struct OtherQuantityHandler::eqn_terms *oqty_terms){
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    const len_t id_Sp = eqsys->GetUnknownID(OptionConstants::UQTY_S_PARTICLE);
    const len_t id_ni = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    const len_t id_nre = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    const len_t id_fhot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    std::string desc = "S_particle = dn_free/dt";

    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));

    // FREE ELECTRON TERM
    Op2->AddTerm(new FreeElectronDensityTransientTerm(fluidGrid,eqsys->GetIonHandler(),id_ni));    
    eqsys->SetOperator(id_Sp, id_ni, Op2);

    // N_RE SOURCES
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
    OptionConstants::eqterm_avalanche_mode ava_mode = 
        (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) {
        // Kinetic 
        real_t pCutoff = s->GetReal("eqsys/n_re/pCutAvalanche");
        Op3->AddTerm(new TotalElectronDensityFromKineticAvalanche(
            fluidGrid, pCutoff, eqsys->GetUnknownHandler(), -1.0
        ));
        desc += " - Gamma_ava*n_re";
    } else if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID || ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_FLUID_HESSLOW) {
        Op3->AddTerm(new AvalancheGrowthTerm(
            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0
        ));
        desc += " - Gamma_ava*n_re";
    }
    OptionConstants::eqterm_compton_mode compton_mode = (enum OptionConstants::eqterm_compton_mode)s->GetInteger("eqsys/n_re/compton/mode");
    if (compton_mode == OptionConstants::EQTERM_COMPTON_MODE_FLUID){
        Op3->AddTerm(new ComptonRateTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );
        desc += " - compton";
    }

    bool hasNreTransport = ConstructTransportTerm(
        Op3, "eqsys/n_re", fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, 
        eqsys->GetUnknownHandler(),s, false, false,
        &oqty_terms->n_re_advective_bc, &oqty_terms->n_re_diffusive_bc
    );
    if(hasNreTransport)
        desc += " - re transport";

    eqsys->SetOperator(id_Sp, id_nre, Op3);

    // F_HOT TRANSPORT TERM
    FVM::Operator *Op_fhot = new FVM::Operator(eqsys->GetHotTailGrid()); // add all kinetic terms not conserving local electron density in this operator
    // Add transport term
    bool hasFHotTransport = ConstructTransportTerm(
        Op_fhot, "eqsys/f_hot", eqsys->GetHotTailGrid(),
        OptionConstants::MOMENTUMGRID_TYPE_PXI, 
        eqsys->GetUnknownHandler(),s, true, false,
        &oqty_terms->f_hot_advective_bc, &oqty_terms->f_hot_diffusive_bc
    );
    if(hasFHotTransport){
        FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
           Op4->AddTerm(new KineticEquationTermIntegratedOverMomentum(
               fluidGrid, eqsys->GetHotTailGrid(), Op_fhot, id_fhot, eqsys->GetUnknownHandler()
           ));
        desc += " - f_hot transport";
        eqsys->SetOperator(id_Sp, id_fhot, Op4);
    }

    eqsys->SetOperator(id_Sp, id_Sp, Op1, desc);
}




