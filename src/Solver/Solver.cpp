/**
 * Implementation of common routines for the 'Solver' routines.
 */

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
Solver::Solver(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations
)
    : unknowns(unknowns), unknown_equations(unknown_equations) {
}


/**
 * Build a jacobian matrix for the equation system.
 *
 * t:    Time to build the jacobian matrix for.
 * mat:  Matrix to use for storing the jacobian.
 * rhs:  Right-hand-side in equation.
 */
void Solver::BuildMatrix(const real_t t, const real_t dt, FVM::BlockMatrix *mat, real_t *S) {
    const len_t nunknowns = unknowns->Size();

    // Build matrix
    mat->Zero();
    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        UnknownQuantityEquation *eqn = unknown_equations->at(i);

        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            mat->SelectSubEquation(nontrivial_unknowns[i], it->first);
            it->second->SetMatrixElements(mat, S);
        }
    }
}

/**
 * Initialize this solver.
 *
 * size:     Number of elements in full unknown vector.
 *           (==> jacobian is of size 'size-by-size').
 * unknowns: List of indices of unknowns to include in the
 *           function vectors/matrices.
 */
void Solver::Initialize(const len_t size, vector<len_t>& unknowns) {
    // Copy list of non-trivial unknowns (those which will
    // appear in the matrices that are built later on)
    nontrivial_unknowns = unknowns;
}

/**
 * Rebuild all equation terms in the equation system for
 * the specified time.
 *
 * t:  Time for which to rebuild the equation system.
 * dt: Length of time step to take next.
 */
void Solver::RebuildTerms(const real_t t, const real_t dt) {
    // Update prescribed quantities
    const len_t N = unknowns->Size();
    for (len_t i = 0; i < N; i++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(i);
        UnknownQuantityEquation *eqn = unknown_equations->at(i);

        if (eqn->IsPrescribed()) {
            eqn->RebuildEquations(t, dt, unknowns);
            FVM::PrescribedParameter *pp = eqn->GetPrescribed();
            uqty->Store(pp->GetData());
        }
    }

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        UnknownQuantityEquation *eqn = unknown_equations->at(i);

        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            it->second->RebuildTerms(t, dt, unknowns);
        }
    }
}

