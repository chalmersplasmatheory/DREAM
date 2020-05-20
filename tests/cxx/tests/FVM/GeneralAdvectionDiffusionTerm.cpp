/**
 * Implementation of a general combined advection-diffusion term
 * for testing purposes.
 */

#include "FVM/Equation/AdvectionDiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "GeneralAdvectionTerm.hpp"
#include "GeneralDiffusionTerm.hpp"
#include "GeneralAdvectionDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;


/**
 * Constructor.
 */
GeneralAdvectionDiffusionTerm::GeneralAdvectionDiffusionTerm(
    DREAM::FVM::Grid *g
) : DREAM::FVM::AdvectionDiffusionTerm(g) {
    
    this->Add(new GeneralAdvectionTerm(g, 1.0));
    this->Add(new GeneralDiffusionTerm(g, 1.0));
}

