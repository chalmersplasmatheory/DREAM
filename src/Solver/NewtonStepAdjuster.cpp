/**
 * Base class for Newton solver step adjustment.
 */

#include "DREAM/Solver/NewtonStepAdjuster.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
NewtonStepAdjuster::NewtonStepAdjuster(
	std::vector<len_t> &nu, FVM::UnknownQuantityHandler *uqh,
	const len_t vs
) : nontrivial_unknowns(nu), unknowns(uqh), vector_size(vs) {
	
	this->x0 = new real_t[vs];
}


/**
 * Destructor.
 */
NewtonStepAdjuster::~NewtonStepAdjuster() {
	delete [] this->x0;
}


/**
 * This routine is called immediately after the Newton step 'dx' has
 * been calculated, and is used to evaluate the new solution 'x'.
 */
void NewtonStepAdjuster::SetX0(
	const len_t, const real_t *x0, const real_t *dx
) {
	for (len_t i = 0; i < this->vector_size; i++)
		this->x0[i] = x0[i];
	
	this->dx = dx;
}


/**
 * Helper routine for updating a Newton solution with a
 * specified damping factor.
 */
void NewtonStepAdjuster::UpdateSolution(real_t *sol, const real_t damping) {
	for (len_t i = 0; i < this->vector_size; i++)
		sol[i] = this->x0[i] - damping * this->dx[i];
}


