/**
 * Definition of equations relating to the hot electron
 * distribution function.
 */

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
        //eqn->AddBoundaryCondition(new FVM::BC::XiInternalBoundaryCondition(hottailGrid));
        // TODO replace this condition with a source term
        //eqn->AddBoundaryCondition(new FVM::BC::PInternalBoundaryCondition(hottailGrid));
    }

    eqsys->SetEquation(OptionConstants::UQTY_F_HOT, OptionConstants::UQTY_F_HOT, eqn);

    // Set initial value of 'f_hot'
    FVM::Interpolator3D *interp = LoadDataR2P(MODULENAME, s, "init");
    enum FVM::Interpolator3D::momentumgrid_type momtype = GetInterp3DMomentumGridType(eqsys->GetHotTailGridType());
    const real_t *init = interp->Eval(hottailGrid, momtype);

    eqsys->SetInitialValue(OptionConstants::UQTY_F_HOT, init);

    delete [] init;
    delete interp;
}

