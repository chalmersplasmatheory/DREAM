/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/n_re"


/**
 * Construct the equation for the runaway electron density, 'n_re'.
 * If the runaway grid is disabled, then n_re is calculated from the
 * particle flux escaping the hot-tail grid, plus any other runaway
 * sources that are enabled.
 */
void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);

    // Add flux from hot tail grid
    if (hottailGrid) {
        // Add the transient term
        FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
        Op_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));

        // Add avalanche growth rate
        Op_nRE->AddTerm(new AvalancheGrowthTerm(fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid(),-1.0) );
        eqsys->SetOperator(id_n_re, id_n_re, Op_nRE);
    

        FVM::Operator *Op_nRE_fHot = new FVM::Operator(fluidGrid);
        len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

        if (eqsys->GetHotTailGridType() == OptionConstants::MOMENTUMGRID_TYPE_PXI) {
            // NOTE We assume that the flux appearing in the equation for 'f_hot'
            // only appears in the (f_hot, f_hot) part of the equation, i.e. in
            // the diagonal block.
            const FVM::Operator *Op = eqsys->GetEquation(id_f_hot)->GetEquation(id_f_hot);

            enum FVM::BC::PXiExternalLoss::bc_type bc =
                (enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger("eqsys/f_hot/boundarycondition");
            Op_nRE_fHot->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
                fluidGrid, Op, id_f_hot, id_n_re, hottailGrid,
                FVM::BC::PXiExternalLoss::BOUNDARY_FLUID, bc
            ));
        } else
            throw NotImplementedException(
                "Currently, the 'DensityFromBoundaryFlux' term only supports "
                "p/xi momentum grids. Hence, you should use a p/xi grid for the "
                "hot-tail distribution function."
            );

        eqsys->SetOperator(id_n_re, id_f_hot, Op_nRE_fHot, "n_re = [flux from f_hot] + n_re*Gamma_ava");


        // Initialize to zero
        eqsys->SetInitialValue(id_n_re, nullptr);
    } else {
        FVM::Operator *Op_nRE = new FVM::Operator(fluidGrid);
        Op_nRE->AddTerm(new FVM::ConstantParameter(fluidGrid, 0));
        eqsys->SetOperator(id_n_re,id_n_re,Op_nRE, "zero");
        eqsys->initializer->AddRule(
            id_n_re,
            EqsysInitializer::INITRULE_EVAL_EQUATION
        );

    }
    
}

