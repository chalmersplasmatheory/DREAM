/**
 * Interface to the Scalable Non-linear Equation Solver (SNES, part of PETSc).
 */

#include <vector>
#include "DREAM/Solver/SolverSNES.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
SolverSNES::SolverSNES(
    FVM::UnknownQuantityHandler *unknowns, 
    vector<UnknownQuantityEquation*> *unknonw_equations
) : Solver(unknowns, unknown_equations) {
}

/**
 * Destructor.
 */
SolverSNES::~SolverSNES() {
    if (F != nullptr)
        delete [] F;
    if (jacobian != nullptr)
        delete jacobian;
}

/**
 * Initialize this solver.
 *
 * size:                Number of elements in full unknown vector.
 *                      (==> jacobian is of size 'size-by-size').
 * nontrivial_unknowns: List of indices of unknowns to include in the
 *                      function vectors/matrices.
 */
void SolverSNES::Initialize(const len_t size, vector<len_t>& nontrivial_unknowns) {
    this->Solver::Initialize(size, nontrivial_unknowns);

    F = new real_t[size];
    jacobian = new FVM::BlockMatrix();

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        UnknownQuantityEquation *eqn = this->unknown_equations->at(i);

        jacobian->CreateSubEquation(eqn->NumberOfElements(), eqn->NumberOfNonZeros_jac());
    }

    jacobian->ConstructSystem();
}

/**
 * Solve the given stage for the non-linear
 * equation system.
 *
 * stage: Stage of the equation system to solve.
 */
void SolverSNES::Solve(const real_t t, const real_t dt) {
    // TODO
    /*RebuildTerms(t, dt);
    BuildVector(F);
    BuildJacobian(jacobian);*/
}

