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
    s->DefineSetting(MODULENAME "/boundarycondition", "Type of boundary condition to use when f_RE is disabled.", (int_t)FVM::BC::PXiExternalLoss::BC_PHI_CONST);
    s->DefineSetting(MODULENAME "/adv_interp/r", "Type of interpolation method to use in r-component of advection term of f_hot kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(MODULENAME "/adv_interp/p1", "Type of interpolation method to use in p1-component of advection term of f_hot kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(MODULENAME "/adv_interp/p2", "Type of interpolation method to use in p2-component of advection term of f_hot kinetic equation.", (int_t)FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED);
    s->DefineSetting(MODULENAME "/adv_interp/fluxlimiterdamping", "Underrelaxation parameter that may be needed to achieve convergence with flux limiter methods", (real_t) 1.0);
    s->DefineSetting(MODULENAME "/pThreshold", "Threshold momentum that defines n_hot from f_hot when resolving thermal population on grid.", (real_t) 10.0);
    s->DefineSetting(MODULENAME "/pThresholdMode", "Unit of provided threshold momentum pThreshold (thermal or mc).", (int_t) FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_THERMAL);
    DefineDataR(MODULENAME,   s, "n0");
    DefineDataR(MODULENAME,   s, "T0");
    DefineDataR2P(MODULENAME, s, "init");
    DefineOptions_Transport(MODULENAME, s, true);
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
    FVM::Operator *eqn = new FVM::Operator(hottailGrid);

    // Add transient term
    eqn->AddTerm(new FVM::TransientTerm(hottailGrid, id_f_hot) );

    string desc;
    // Determine whether electric field acceleration should be
    // modelled with an advection or a diffusion term
    //
    // XXX Here we assume that all momentum grids have
    // the same grid points
    if (eqsys->GetHotTailGridType() == OptionConstants::MOMENTUMGRID_TYPE_PXI &&
            hottailGrid->GetMomentumGrid(0)->GetNp2() == 1) {
        
        desc = "Reduced kinetic equation";

        eqn->AddTerm(new ElectricFieldDiffusionTerm(
            hottailGrid, eqsys->GetHotTailCollisionHandler(), eqsys->GetUnknownHandler())
        );
    // Model as an advection term
    } else {
        desc = "3D kinetic equation";

        // Electric field term
        eqn->AddTerm(new ElectricFieldTerm(
            hottailGrid, eqsys->GetUnknownHandler(), eqsys->GetHotTailGridType()
        ));

        // Pitch scattering term
        eqn->AddTerm(new PitchScatterTerm(
            hottailGrid, eqsys->GetHotTailCollisionHandler(), eqsys->GetHotTailGridType(),
            eqsys->GetUnknownHandler()
        ));
    }

    // ALWAYS PRESENT

    // Slowing down term
    eqn->AddTerm(new SlowingDownTerm(
            hottailGrid, eqsys->GetHotTailCollisionHandler(), eqsys->GetHotTailGridType(), 
            eqsys->GetUnknownHandler()
        )
    );

    // Energy diffusion
    eqn->AddTerm(new EnergyDiffusionTerm(
        hottailGrid, eqsys->GetHotTailCollisionHandler(), eqsys->GetHotTailGridType(),
        eqsys->GetUnknownHandler()
    ));
    
    // Add transport term
    ConstructTransportTerm(
        eqn, MODULENAME, hottailGrid,
        eqsys->GetHotTailGridType(), s, true
    );

    // EXTERNAL BOUNDARY CONDITIONS
    // Lose particles to runaway region
	if (eqsys->HasRunawayGrid()) {
		len_t id_f_re = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);
		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalKineticKinetic(
			hottailGrid, hottailGrid, eqsys->GetRunawayGrid(), eqn,
			id_f_hot, id_f_re, FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
		));
	} else {
		enum FVM::BC::PXiExternalLoss::bc_type bc =
			(enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger(MODULENAME "/boundarycondition");

		eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
			hottailGrid, eqn, id_f_hot, id_f_hot, nullptr,
			FVM::BC::PXiExternalLoss::BOUNDARY_KINETIC, bc
		));
	}

    // Set interpolation scheme for advection term
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_r =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(MODULENAME "/adv_interp/r");
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_p1 =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(MODULENAME "/adv_interp/p1");
    enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_p2 =
			(enum FVM::AdvectionInterpolationCoefficient::adv_interpolation)s->GetInteger(MODULENAME "/adv_interp/p2");
    real_t fluxLimiterDamping = (real_t)s->GetReal(MODULENAME "/adv_interp/fluxlimiterdamping");
    eqn->SetAdvectionInterpolationMethod(adv_interp_r,  FVM::FLUXGRIDTYPE_RADIAL, id_f_hot, fluxLimiterDamping);
    eqn->SetAdvectionInterpolationMethod(adv_interp_p1, FVM::FLUXGRIDTYPE_P1,     id_f_hot, fluxLimiterDamping);
    eqn->SetAdvectionInterpolationMethod(adv_interp_p2, FVM::FLUXGRIDTYPE_P2,     id_f_hot, fluxLimiterDamping);

    eqsys->SetOperator(id_f_hot, id_f_hot, eqn, desc);


    // AVALANCHE SOURCE TERM
    OptionConstants::eqterm_avalanche_mode ava_mode = (enum OptionConstants::eqterm_avalanche_mode)s->GetInteger("eqsys/n_re/avalanche");
    if(ava_mode == OptionConstants::EQTERM_AVALANCHE_MODE_KINETIC){
        len_t id_n_re = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);

        if(eqsys->GetHotTailGridType() != OptionConstants::MOMENTUMGRID_TYPE_PXI)
            throw NotImplementedException("f_hot: Kinetic avalanche source only implemented for p-xi grid.");
    
        real_t pCutoff = s->GetReal("eqsys/n_re/pCutAvalanche");
        FVM::Operator *Op_ava = new FVM::Operator(hottailGrid);
        Op_ava->AddTerm(new AvalancheSourceRP(hottailGrid, eqsys->GetUnknownHandler(), pCutoff, -1.0 ));
        eqsys->SetOperator(id_f_hot, id_n_re, Op_ava);
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
        Op4->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0e-5));
        eqsys->SetOperator(id_Sp, id_n_cold, Op1, "integral(f_hot) = n_cold + n_hot + 1e-5*S_p");
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


    // Set initial value of 'f_hot'
    //   First, we check whether the distribution has been specified numerically.
    //   If it hasn't, we prescribe a Maxwellian with the correct temperature.
    len_t nx[3];
    if (s->GetRealArray(MODULENAME "/init/x", 3, nx, false) != nullptr) {
        FVM::Interpolator3D *interp = LoadDataR2P(MODULENAME, s, "init");
        enum FVM::Interpolator3D::momentumgrid_type momtype = GetInterp3DMomentumGridType(eqsys->GetHotTailGridType());
        const real_t *init = interp->Eval(hottailGrid, momtype);

        eqsys->SetInitialValue(id_f_hot, init);

        delete [] init;
        delete interp;
    } else {
        real_t *n0 = LoadDataR(MODULENAME, hottailGrid->GetRadialGrid(), s, "n0");
        real_t *T0 = LoadDataR(MODULENAME, hottailGrid->GetRadialGrid(), s, "T0");

        ConstructEquation_f_maxwellian(OptionConstants::UQTY_F_HOT, eqsys, hottailGrid, n0, T0);

        delete [] T0;
        delete [] n0;
    }
}

/**
 * Initializes the hot electron distribution function as a
 * Maxwellian with the specified temperature and density
 * profiles.
 *
 * n0: Initial density profile of hot electrons.
 * T0: Initial temperature profile of hot electrons.
 */
void SimulationGenerator::ConstructEquation_f_maxwellian(
    const std::string& uqtyName,
    EquationSystem *eqsys, FVM::Grid *grid, const real_t *n0, const real_t *T0
) {
    const len_t nr = grid->GetNr();
    real_t *init = new real_t[grid->GetNCells()];

    // Construct Maxwellian
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np1 = grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = grid->GetMomentumGrid(ir)->GetNp2();
        const real_t *pvec = grid->GetMomentumGrid(ir)->GetP();

        // Define distribution offset vector
        real_t *f = init + offset;
        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1; i++)
                f[j*np1 + i] = Constants::RelativisticMaxwellian(pvec[j*np1+i], n0[ir], T0[ir]);

        offset += np1*np2;
    }

    eqsys->SetInitialValue(uqtyName, init);

    delete [] init;
}

