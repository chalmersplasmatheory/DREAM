/**
 * Implementation of a "constant" equation parameter.
 */

#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 *
 * g: Grid on which the parameter lives.
 * v: Value of parameter.
 */
ConstantParameter::ConstantParameter(Grid *g, const real_t v)
    : PredeterminedParameter(g) {

    const len_t N = g->GetNCells();

    for (len_t i = 0; i < N; i++)
        this->currentData[i] = v;
}

/**
 * Destructor.
 */
ConstantParameter::~ConstantParameter() { }

