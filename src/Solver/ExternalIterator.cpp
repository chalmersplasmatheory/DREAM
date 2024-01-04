/**
 * Module for further iterating on an equation system. This module is
 * responsible for updating some of the unknowns in the equation system.
 */

#include <vector>
#include "DREAM/Solver/ExternalIterator.hpp"
#include "DREAM/Solver/Solver.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
ExternalIterator::ExternalIterator(
	FVM::UnknownQuantityHandler *unknowns,
	vector<UnknownQuantityEquation*> *unknown_equations,
	const bool verbose
) : unknowns(unknowns), unknown_equations(unknown_equations),
	printVerbose(verbose) { }

/**
 * Destructor.
 */
ExternalIterator::~ExternalIterator() {
	delete [] this->vec;
	delete [] this->dvec;
}


/**
 * Initialize this module.
 *
 * unknowns: List of IDs of the unknowns to be evaluated by this module.
 */
void ExternalIterator::Initialize(
	vector<len_t> &ext_unknowns
) {
	this->external_unknowns = ext_unknowns;
	this->vecsize = this->unknowns->GetLongVectorSize(ext_unknowns);

	this->vec = new real_t[this->vecsize];
	this->dvec = new real_t[this->vecsize];
}

/**
 * Evaluate the quantities handled by this module. Returns true if all external
 * quantities are converged.
 */
bool ExternalIterator::Solve(
	const real_t t, const real_t dt
) {
	len_t vecoffs = 0;
	for (len_t uqnId : external_unknowns) {
		UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);

		if (eqn->GetOperators().size() != 1) {
			throw SolverException(
				"External quantity '%s': Equation must contain exactly one "
				"operator.",
				unknowns->GetUnknown(uqnId)->GetName().c_str()
			);
		}

		// Fetch the operator of the equation
		FVM::Operator *op = nullptr;
		for (auto o : eqn->GetOperators())
			op = o.second;

		if (op == nullptr)
			continue;

		// Rebuild
		op->RebuildTerms(t, dt, unknowns);
		
		const real_t *data = unknowns->GetUnknownData(uqnId);
		op->SetVectorElements(this->vec + vecoffs, data);

		vecoffs += unknowns->GetUnknown(uqnId)->NumberOfElements();
	}

	// Evaluate norm dx = x_new - x_old
	this->unknowns->GetLongVector(this->external_unknowns, this->dvec);
	for (len_t i = 0; i < vecoffs; i++)
		dvec[i] = this->vec[i] - dvec[i];

	// Check if ||dx||^2 is small enough for convergence
	bool conv = this->convChecker->IsConverged(vec, dvec, this->printVerbose);

	// Update solution value
	this->unknowns->Store(this->external_unknowns, this->vec);

	return conv;
}

/**
 * Set the convergence checker to use.
 */
void ExternalIterator::SetConvergenceChecker(
	ConvergenceChecker *cc
) {
	this->convChecker = cc;
}

