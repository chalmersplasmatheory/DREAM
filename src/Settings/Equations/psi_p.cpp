/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

#include <iostream>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"
#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include "FVM/Equation/WeightedTransientTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/AmperesLawDiffusionTerm.hpp"

using namespace DREAM;
using namespace std;

#define MODULENAME "eqsys/psi_p"


/**
 * Define settings for the hot-tail distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_psi_p(Settings *s) {
    DefineDataR2P(MODULENAME, s, "init");
}

/**
 * Construct the equation for the hot electron distribution
 * function. This function is only called if the hot-tail grid
 * is enabled.
 *
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_psi_p(
    EquationSystem *eqsys, Settings * /*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Equation *eqn_E1 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn_E2 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn_j1 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn_j2 = new FVM::Equation(fluidGrid);


    // weightFunc1 represents R0*<B*nabla phi>/Bmin
    std::function<real_t(len_t,len_t,len_t)> weightFunc1 = [fluidGrid](len_t ir,len_t, len_t)
        {return fluidGrid->GetRadialGrid()->GetFSA_1OverR2(ir) * fluidGrid->GetRadialGrid()->GetBTorG(ir) / fluidGrid->GetRadialGrid()->GetBmin(ir);};
    // Add transient term
    eqn_E1->AddTerm(new FVM::WeightedTransientTerm(
        fluidGrid, eqsys->GetUnknownID(OptionConstants::UQTY_POL_FLUX), &weightFunc1)
    );
    // eqn_E1->AddBoundaryCondition(new PsiFluxWall(wallSettings));
    
    
    std::function<real_t(len_t,len_t,len_t)> weightFunc2 = [fluidGrid](len_t ir,len_t, len_t)
        {return 2*M_PI*sqrt(fluidGrid->GetRadialGrid()->GetFSA_B2(ir));};
    eqn_E2->AddTerm(new FVM::WeightedIdentityTerm(
        fluidGrid,&weightFunc2) 
    );


    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_POL_FLUX, eqn_E1, "Poloidal flux resistive diffusion equation");
    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn_E2);
    
    bool settingHyperresistivity = false;
    if(settingHyperresistivity){
        FVM::Equation *eqn_E3 = new FVM::Equation(fluidGrid);
//        eqn_E3->AddTerm(new HyperresistiveTerm(
//           fluidGrid, eqsys->GetTransportHandler(),eqsys->GetUnknownHandler()) 
//       );
        eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, eqn_E3);
    } 


    

    // weightFunc3 represents mu0*R0*<B*nabla phi>/Bmin (ie mu0*weightFunc1)
    std::function<real_t(len_t,len_t,len_t)> weightFunc3 = [fluidGrid](len_t ir,len_t, len_t)
        {return Constants::mu0 * fluidGrid->GetRadialGrid()->GetFSA_1OverR2(ir) * fluidGrid->GetRadialGrid()->GetBTorG(ir) / fluidGrid->GetRadialGrid()->GetBmin(ir);};

    eqn_j1->AddTerm(new FVM::WeightedIdentityTerm(
        fluidGrid, &weightFunc3)
    );

    eqn_j2->AddTerm(new AmperesLawDiffusionTerm(
        fluidGrid)
    );

    /**
     * TODO: Add additional boundary conditions.
     */


    eqsys->SetEquation(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_J_TOT, eqn_j1, "Poloidal flux Ampere's law");
    eqsys->SetEquation(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_POL_FLUX, eqn_j2);

/*
    // Set initial value of 'f_hot'
    //   First, we check whether the distribution has been specified numerically.
    //   If it hasn't, we prescribe a Maxwellian with the correct temperature.
    len_t nx[3];
    if (s->GetRealArray(MODULENAME "/init/x", 3, nx, false) != nullptr) {
        FVM::Interpolator3D *interp = LoadDataR2P(MODULENAME, s, "init");
        enum FVM::Interpolator3D::momentumgrid_type momtype = GetInterp3DMomentumGridType(eqsys->GetHotTailGridType());
        const real_t *init = interp->Eval(hottailGrid, momtype);

        eqsys->SetInitialValue(OptionConstants::UQTY_F_HOT, init);

        delete [] init;
        delete interp;
    } else {
        real_t *n0 = LoadDataR(MODULENAME, hottailGrid->GetRadialGrid(), s, "n0");
        real_t *T0 = LoadDataR(MODULENAME, hottailGrid->GetRadialGrid(), s, "T0");

        ConstructEquation_f_hot_maxwellian(eqsys, hottailGrid, n0, T0);

        delete [] T0;
        delete [] n0;
    }
*/
}
