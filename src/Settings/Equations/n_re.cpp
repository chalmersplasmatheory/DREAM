/**
 * Definition of equations relating to n_re (the radial density
 * of runaway electrons).
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "DREAM/Equations/Fluid/DensityFromDistributionFunction.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "equationsystem/n_re"


/**
 * Construct the equation for the runaway electron density, 'n_re'.
 * If the runaway grid is disabled, then n_re is calculated from the
 * particle flux escaping the hot-tail grid, plus any other runaway
 * sources that are enabled. If the runaway grid is enabled, however,
 * then 'n_re' is calculated by taking the density moment of the
 * runaway electron distribution function, 'f_re'.
 */
void SimulationGenerator::ConstructEquation_n_re(
    EquationSystem *eqsys, Settings*
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();

    len_t id_n_re  = eqsys->GetUnknownID(OptionConstants::UQTY_N_RE);

    // If the runaway grid is enabled, calculate as a moment of 'f_re'...
    if (runawayGrid) {
        len_t id_f_re  = eqsys->GetUnknownID(OptionConstants::UQTY_F_RE);

        FVM::Equation *eqn = new FVM::Equation(fluidGrid);

        DensityFromDistributionFunction *mq = new DensityFromDistributionFunction(
            fluidGrid, runawayGrid, id_n_re, id_f_re
        );
        eqn->AddTerm(mq);
        eqsys->SetEquation(id_n_re, id_f_re, eqn, "Moment of f_re");
    // Otherwise, as flux of particles from the hot-tail grid + source terms
    } else {
        if (hottailGrid) {
            FVM::Equation *eqn_nRE_fHot = new FVM::Equation(fluidGrid);
            len_t id_f_hot = eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT);

            if (eqsys->GetHotTailGridType() == OptionConstants::MOMENTUMGRID_TYPE_PXI) {
                // NOTE We assume that the flux appearing in the equation for 'f_hot'
                // only appears in the (f_hot, f_hot) part of the equation, i.e. in
                // the diagonal block.
                const FVM::Equation *eqn = eqsys->GetEquation(id_f_hot)->GetEquation(id_f_hot);
                DensityFromBoundaryFluxPXI *mq = new DensityFromBoundaryFluxPXI(
                    fluidGrid, hottailGrid, eqn
                );

                eqn_nRE_fHot->AddTerm(mq);
            } else
                throw NotImplementedException(
                    "Currently, the 'DensityFromBoundaryFlux' term only supports "
                    "p/xi momentum grids. Hence, you should use a p/xi grid for the "
                    "hot-tail distribution function."
                );

            eqsys->SetEquation(id_n_re, id_f_hot, eqn_nRE_fHot, "Flux from f_hot");
        }

        // TODO Other source terms

        // Add the transient term
        FVM::Equation *eqn_nRE = new FVM::Equation(fluidGrid);
        eqn_nRE->AddTerm(new FVM::TransientTerm(fluidGrid, id_n_re));

        eqsys->SetEquation(id_n_re, id_n_re, eqn_nRE);
        
        // Initialize to zero
        eqsys->SetInitialValue(id_n_re, nullptr);
    }
}

