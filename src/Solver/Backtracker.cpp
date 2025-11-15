/**
 * This file contains an implementation of backtracking for the Newton solver.
 * Backtracking can be used to accelerate the Newton solver when starting from
 * a solution that is far from the true solution.
 *
 * Notes on the implementation of this module can be found in
 * 'doc/notes/discretisation.pdf'.
 */

#include <string>
#include <vector>
#include <petsc.h>
#include "DREAM/Solver/Backtracker.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
Backtracker::Backtracker(
    std::vector<len_t>& nu, FVM::UnknownQuantityHandler *uqh,
    IonHandler *ions, std::vector<std::string>& monitor
) : PhysicalStepAdjuster(nu, uqh, ions) {

    this->nnu = nu.size();

	// Set up list of unknowns to monitor
	this->nMonitor = monitor.size();
	if (this->nMonitor > 0) {
		this->monitor = new len_t[this->nMonitor];
		for (len_t i = 0; i < this->nMonitor; i++)
			this->monitor[i] = uqh->GetUnknownID(monitor[i]);
	} else {
		// If the list is empty, we consider all non-trivial unknowns
		this->nMonitor = nnu;
		this->monitor = new len_t[this->nMonitor];
		for (len_t i = 0; i < this->nMonitor; i++)
			this->monitor[i] = nu[i];
	}

    this->f0 = new real_t[nnu];
    this->f1 = new real_t[nnu];
    this->f2 = new real_t[nnu];
    this->gradf_deltax = new real_t[nnu];
    this->decreasing = new bool[nnu];
}

/**
 * Destructor.
 */
Backtracker::~Backtracker() {
    delete [] this->decreasing;
    delete [] this->gradf_deltax;
    delete [] this->f2;
    delete [] this->f1;
    delete [] this->f0;
	delete [] this->monitor;
}


/**
 * Determine the factor 'lambda' with which the next Newton step
 * should be adjusted.
 */
real_t Backtracker::Adjust(
    len_t iteration,
    const real_t *x, const real_t *dx, Vec &F,
    FVM::BlockMatrix *jac
) {
    real_t lambda=1;

    EvaluateTargetFunction(F, jac);

    // Is backtracking necessary?
    if (!IsDecreasing(decreasing)) {
        if (this->nIteration > 1)
            lambda = CalculateLambda();
        else {} // First, we take a step with lambda=1
    } else
        // Successfull Newton step -- reset backtracking
        ResetBacktracking();
    
    this->lambda2 = this->lambda1;
    this->lambda1 = lambda;

    // Prevent Newton step from making key quantities negative
    real_t alpha = this->PhysicalStepAdjuster::Adjust(iteration, x, dx, F, jac);

    // If step must be constrained further due to physical reasons,
    // do so...
    if (lambda > alpha)
        lambda = alpha;

    return lambda;
}

/**
 * Calculate the helper function g(\lambda).
 */
real_t Backtracker::CalculateLambda() {
    real_t lambda = 1;

    // Iterate over monitored unknowns
    for (len_t i = 0; i < this->nMonitor; i++) {
        real_t l = 1;

        if (this->nIteration == 2) {    // quadratic approximation to g(lambda)
            real_t gp0 = gradf_deltax[i];
            real_t g1  = f1[i];
            real_t g0  = f0[i];

            l = -gp0 / (2*(g1 - g0 - gp0));
        } else if (this->nIteration > 2) {  // cubic approximation to g(lambda)
            real_t gp0 = gradf_deltax[i];
            real_t g0  = f0[i];
            real_t g1  = f1[i];
            real_t g2  = f2[i];
            real_t l1  = lambda1, l2 = lambda2;

            real_t l1_2 = l1*l1, l1_3 = l1_2*l1;
            real_t l2_2 = l2*l2, l2_3 = l2_2*l2;

            real_t den = l1_2*l2_2*(l1-l2);
            real_t a = (l2_3*g1 - l1_3*g2 + l1*l2*(l1-l2)*gp0+ (l1_2-l2_2)*g0) / den;
            real_t b = (l1_3*g2 - l2_3*g1 + l1*l2*(l2_3-l1_3)*gp0 + (l2_3 - l1_3)*g0) / den;

            l = (-b + sqrt(b*b - 3*a*gp0)) / (3*a);
        }

        // Select the smallest (most limited) lambda for
        // each of the unknowns...
        if (l < lambda) {
            lambda = l;
            this->limitingUnknown = monitor[i];
        }
    }

    // Prevent slow convergence by limiting lambda...
    if (lambda < 0.1*lambda1)
        lambda = 0.1*lambda;
    else if (lambda > 0.5*lambda)
        lambda = 0.5*lambda;
    
    return lambda;
}

/**
 * Evaluate the target function f = -F*F, where F is
 * the residual function.
 */
void Backtracker::EvaluateTargetFunction(Vec &F, FVM::BlockMatrix *jac) {
    real_t *fvec;

    VecGetArray(F, &fvec);

    for (len_t i = 0; i < this->nMonitor; i++) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(monitor[i]);
        len_t offset = jac->GetOffsetById(monitor[i]);
        len_t N = uqn->NumberOfElements();

        // Copy value of f(x) at previous lambda
        f2[i] = f1[i];

        f1[i] = 0;
        for (len_t j = 0; j < N; j++)
            f1[i] += fvec[j+offset]*fvec[j+offset];

        // On the very first backtracking iteration, also calculate
        // the gradient of f, dotted with the Newton step 'delta-x',
        // and store g(0) = f
        if (nIteration == 0) {
            gradf_deltax[i] = -f1[i];
            f0[i] = 0.5*f1[i];
        }

        f1[i] *= 0.5;
    }

    VecRestoreArray(F, &fvec);

    nIteration++;
}

/**
 * Check whether the most recent step satisfies the condition for
 * the step to decrease the target function f. Returns 'true' if
 * the target function for each unknown is decreasing. Otherwise,
 * detailed results are available in 'dec'.
 *
 * dec: Each element corresponds to a non-trivial unknown of the
 *      system. If the element is 'true', the Newton step is leading
 *      to a more optimal value for the unknown.
 */
bool Backtracker::IsDecreasing(bool *dec) {
    // Assume that all is well on the first Newton step...
    if (this->nIteration == 0)
        return true;

    bool v = true;
    for (len_t i = 0; i < this->nMonitor; i++) {
        dec[i] = (f1[i] > f2[i]+ALPHA*gradf_deltax[i]);
        v = v && dec[i];
    }

    return v;
}

/**
 * Reset the backtracking algorithm in preparation for the
 * next Newton iteration.
 */
void Backtracker::ResetBacktracking() {
    this->nIteration = 0;
}

