/**
 * A collection of fluid runaway source terms.
 */

#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
RunawaySourceTerm::RunawaySourceTerm(FVM::Grid *grid)
    : EquationTerm(grid) {
}

