/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

#include <gsl/gsl_sf_bessel.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldTerm.hpp"
#include "DREAM/Equations/Kinetic/ElectricFieldDiffusionTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Equation/BoundaryConditions/PInternalBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryConditions/XiInternalBoundaryCondition.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Interpolator3D.hpp"


using namespace DREAM;

#define MODULENAME "equationsystem/f_hot"


/**
 * Define settings for the hot-tail distribution function.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_f_hot(Settings *s) {
    DefineDataR(MODULENAME, s, "n0");
    DefineDataR(MODULENAME, s, "T0");
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
void SimulationGenerator::ConstructEquation_f_hot(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::Equation *eqn = new FVM::Equation(hottailGrid);

    // Add transient term
    FVM::TransientTerm *tt = new FVM::TransientTerm(hottailGrid, eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT));
    eqn->AddTerm(tt);

    // Determine whether electric field acceleration should be
    // modelled with an advection or a diffusion term
    //
    // XXX Here we assume that all momentum grids have
    // the same grid points
    if (eqsys->GetHotTailGridType() == OptionConstants::MOMENTUMGRID_TYPE_PXI &&
        hottailGrid->GetMomentumGrid(0)->GetNp2() == 1) {
        
        throw SettingsException(
            "f_hot: No support implemented for collisions yet, so E-field cannot be modelled "
            "as a diffusion term. Please set nXi > 1."
        );

        ElectricFieldDiffusionTerm *efdt = new ElectricFieldDiffusionTerm(
            hottailGrid, nullptr, eqsys->GetUnknownHandler()
        );

        eqn->AddTerm(efdt);

        // BOUNDARY CONDITIONS
        // Lose particles to runaway region
        eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(hottailGrid, eqn));
        // Standard internal boundary conditions
        eqn->AddBoundaryCondition(new FVM::BC::XiInternalBoundaryCondition(hottailGrid));
        // TODO replace this condition with a source term
        eqn->AddBoundaryCondition(new FVM::BC::PInternalBoundaryCondition(hottailGrid));
        
    // Model as an advection term
    } else {
        ElectricFieldTerm *eft = new ElectricFieldTerm(hottailGrid, eqsys->GetUnknownHandler(), eqsys->GetHotTailGridType());
        eqn->AddTerm(eft);

        // BOUNDARY CONDITIONS
        // Lose particles to runaway region
        eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(hottailGrid, eqn));
        // Standard internal boundary conditions
        eqn->AddBoundaryCondition(new FVM::BC::XiInternalBoundaryCondition(hottailGrid));
        // TODO replace this condition with a source term
        eqn->AddBoundaryCondition(new FVM::BC::PInternalBoundaryCondition(hottailGrid));
    }

    eqsys->SetEquation(OptionConstants::UQTY_F_HOT, OptionConstants::UQTY_F_HOT, eqn);

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
}

/**
 * Initializes the hot electron distribution function as a
 * Maxwellian with the specified temperature and density
 * profiles.
 *
 * n0: Initial density profile of hot electrons.
 * T0: Initial temperature profile of hot electrons.
 */
void SimulationGenerator::ConstructEquation_f_hot_maxwellian(
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
        // Normalized temperature and scale factor
        real_t Theta  = T0[ir] / Constants::mc2inEV;
        real_t tK2exp = Theta * gsl_sf_bessel_Knu_scaled(2.0, 1.0/Theta);

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                const real_t p = pvec[j*np1+i];
                const real_t g = sqrt(1+p*p);

                f[j*np1 + i] = g*p*n0[ir] / tK2exp * exp((1-g)/Theta);
            }
        }

        offset += np1*np2;
    }

    eqsys->SetInitialValue(OptionConstants::UQTY_F_HOT, init);

    delete [] init;
}

