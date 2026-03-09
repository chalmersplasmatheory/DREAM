/**
 * Implementation of a rescaled electric field term, for use in Ohm's law.
 */

#include <cmath>
#include "DREAM/Equations/Fluid/OhmicElectricFieldTerm.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
OhmicElectricFieldTerm::OhmicElectricFieldTerm(
	FVM::Grid *grid, const real_t scaleFactor
) : DiagonalLinearTerm(grid), scaleFactor(scaleFactor) { }


/**
 * Set the weights of this term.
 */
void OhmicElectricFieldTerm::SetWeights() {
	len_t nr = grid->GetNr();
	FVM::RadialGrid *rgrid = grid->GetRadialGrid();
	for (len_t ir = 0; ir < nr; ir++)
		weights[ir] = scaleFactor * std::sqrt(rgrid->GetFSA_B2(ir));
}

