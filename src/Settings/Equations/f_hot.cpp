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
    s->DefineSetting(MODULENAME "/pThreshold", "Threshold momentum that defines n_hot from f_hot when resolving thermal population on grid.", (real_t) 5.0);
    s->DefineSetting(MODULENAME "/pThresholdMode", "Unit of provided threshold momentum pThreshold (thermal or mc).", (int_t) FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL);
    s->DefineSetting(MODULENAME "/particleSource", "Include particle source which enforces the integral over the distribution to follow n_hot+n_cold.", (int_t) OptionConstants::EQTERM_PARTICLE_SOURCE_EXPLICIT);
    s->DefineSetting(MODULENAME "/particleSourceShape", "Determines the shape of the particle source term.", (int_t)OptionConstants::EQTERM_PARTICLE_SOURCE_SHAPE_MAXWELLIAN);

    s->DefineSetting(MODULENAME "/dist_mode", "Which analytic model to use for the hottail distribution", (int_t)OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL);
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
        &oqty_terms->f_hot_advective_bc, &oqty_terms->f_hot_diffusive_bc,
        &oqty_terms->f_hot_ripple_Dxx, rescaleMaxwellian
    );

    // Add kinetic-kinetic boundary condition if necessary...
    if (!addExternalBC) {
		len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			hottailGrid, hottailGrid, eqsys->GetRunawayGrid(), eqn,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
		));
    }

    // Add avalanche source
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) {
        if(eqsys->GetHotTailGridType() != OptionConstants::MOMENTUMGRID_TYPE_PXI)
            throw NotImplementedException("f_hot: Kinetic avalanche source only implemented for p-xi grid.");

        real_t pCutoff = s->GetReal("eqsys/n_re/pCutAvalanche");
        FVM::Operator *Op_ava = new FVM::Operator(hottailGrid);
        Op_ava->AddTerm(new AvalancheSourceRP(hottailGrid, eqsys->GetUnknownHandler(), pCutoff, -1.0 ));
        len_t id_n_re = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_N_RE);
        eqsys->SetOperator(id_f_hot, id_n_re, Op_ava);
    }

    // PARTICLE SOURCE TERMS
    const len_t id_Sp = eqsys->GetUnknownID(OptionConstants::UQTY_S_PARTICLE);
    FVM::Operator *Op_source = new FVM::Operator(hottailGrid);

    // Enable particle source term ?
    OptionConstants::eqterm_particle_source_mode particleSource = (OptionConstants::eqterm_particle_source_mode)s->GetInteger(MODULENAME "/particleSource"); 
    OptionConstants::eqterm_particle_source_shape pSourceShape  = (OptionConstants::eqterm_particle_source_shape)s->GetInteger(MODULENAME "/particleSourceShape");

    enum ParticleSourceTerm::ParticleSourceShape shape;
    switch (pSourceShape) {
        case OptionConstants::EQTERM_PARTICLE_SOURCE_SHAPE_MAXWELLIAN:
            shape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_MAXWELLIAN;
            break;

        case OptionConstants::EQTERM_PARTICLE_SOURCE_SHAPE_DELTA:
            shape = ParticleSourceTerm::PARTICLE_SOURCE_SHAPE_DELTA;
            break;
        
        default:
            throw SettingsException("Unrecognized particle source term shape: %d.", pSourceShape);
    }

    // Construct particle source term
    Op_source->AddTerm(
        new ParticleSourceTerm(
            hottailGrid,eqsys->GetUnknownHandler(), shape
        )
    );
    eqsys->SetOperator(id_f_hot, id_Sp, Op_source);

    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    OptionConstants::collqty_collfreq_mode collfreq_mode = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
    if(particleSource==OptionConstants::EQTERM_PARTICLE_SOURCE_IMPLICIT && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL))
        ConstructEquation_S_particle_implicit(eqsys, s);
    else if(particleSource==OptionConstants::EQTERM_PARTICLE_SOURCE_EXPLICIT && (collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL))
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
        real_t pLower, pUpper, scaleFactor;
        TotalElectronDensityFromKineticAvalanche(FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, real_t scaleFactor = 1.0) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_TOT),u), pLower(pLower), pUpper(pUpper), scaleFactor(scaleFactor) {}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = scaleFactor * AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(pLower, pUpper);
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
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    
    const len_t id_Sp = eqsys->GetUnknownID(OptionConstants::UQTY_S_PARTICLE);
    const len_t id_ni = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    const len_t id_nre = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);
    const len_t id_ntot = eqsys->GetUnknownID(OptionConstants::UQTY_N_TOT);
    const len_t id_fhot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

    std::string desc = "S_particle = dn_free/dt";

    FVM::Operator *Op_Sp = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Ni = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Ntot = new FVM::Operator(fluidGrid);
    FVM::Operator *Op_Nre = new FVM::Operator(fluidGrid);

    Op_Sp->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0));

    // FREE ELECTRON TERM
    Op_Ni->AddTerm(new FreeElectronDensityTransientTerm(fluidGrid,eqsys->GetIonHandler(),id_ni));    
    eqsys->SetOperator(id_Sp, id_ni, Op_Ni);

    // N_RE SOURCES
    
    // Add contribution from kinetic avalanche source
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC) {
        if(eqsys->GetHotTailGridType() != OptionConstants::MOMENTUMGRID_TYPE_PXI)
            throw NotImplementedException("f_hot: Kinetic avalanche source only implemented for p-xi grid.");

        real_t pCutoff = s->GetReal("eqsys/n_re/pCutAvalanche");
        real_t pMax = hottailGrid->GetMomentumGrid(0)->GetP1_f(hottailGrid->GetNp1(0));
        Op_Nre->AddTerm(
            new TotalElectronDensityFromKineticAvalanche(fluidGrid, pCutoff, pMax, unknowns, -1.0)
        );
        desc += " - internal avalanche";
    }
    // Add source terms
    bool signPositive = false;
    RunawaySourceTermHandler *rsth = ConstructRunawaySourceTermHandler(
        fluidGrid, hottailGrid, eqsys->GetRunawayGrid(), fluidGrid, eqsys->GetUnknownHandler(),
        eqsys->GetREFluid(), eqsys->GetIonHandler(), eqsys->GetAnalyticHottailDistribution(), oqty_terms, s, signPositive
    );

    rsth->AddToOperators(Op_Nre, Op_Ntot, Op_Ni);
    desc += rsth->GetDescription();

    bool hasNreTransport = ConstructTransportTerm(
        Op_Nre, "eqsys/n_re", fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, 
        eqsys->GetUnknownHandler(),s, false, false,
        &oqty_terms->n_re_advective_bc, &oqty_terms->n_re_diffusive_bc
    );
    if(hasNreTransport)
        desc += " - re transport";

    eqsys->SetOperator(id_Sp, id_nre, Op_Nre);
    eqsys->SetOperator(id_Sp, id_ntot, Op_Ntot);

    // F_HOT TRANSPORT TERM
    FVM::Operator *Op_fhot_tmp = new FVM::Operator(eqsys->GetHotTailGrid()); // add all kinetic terms not conserving local electron density in this operator
    // Add transport term
    bool hasFHotTerm = ConstructTransportTerm(
        Op_fhot_tmp, "eqsys/f_hot", eqsys->GetHotTailGrid(),
        OptionConstants::MOMENTUMGRID_TYPE_PXI, 
        eqsys->GetUnknownHandler(),s, true, false,
        &oqty_terms->f_hot_advective_bc, &oqty_terms->f_hot_diffusive_bc
    );
    if(hasFHotTerm){ // add kinetic term integrated over momentum
        FVM::Operator *Op_fhot = new FVM::Operator(fluidGrid);
           Op_fhot->AddTerm(new KineticEquationTermIntegratedOverMomentum(
               fluidGrid, eqsys->GetHotTailGrid(), Op_fhot_tmp, id_fhot, eqsys->GetUnknownHandler()
           ));
        desc += " - f_hot transport";
        eqsys->SetOperator(id_Sp, id_fhot, Op_fhot);
    }

    eqsys->SetOperator(id_Sp, id_Sp, Op_Sp, desc);
}




