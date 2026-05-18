/**
 * A term representing the RHS of
 *   
 *   T_hot = (2/3) * W_hot * n_hot
 */

#include "DREAM/Equations/Fluid/HotElectronTemperatureTerm.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
HotElectronTemperatureTerm::HotElectronTemperatureTerm(
	FVM::Grid *g, FVM::UnknownQuantityHandler *u
) : FVM::EvaluableEquationTerm(g), unknowns(u) {
}


/**
 * Transform to apply to evaluate this equation term.
 */



