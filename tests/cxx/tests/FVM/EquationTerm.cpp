/**
 * Implementation of general tests for all 'EquationTerm' classes.
 */

#include <petsc.h>
#include "FVM/Matrix.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "EquationTerm.hpp"


using namespace DREAMTESTS::FVM;


/**
 * Check wether this equation term is conservative
 * on all the available grids.
 */
bool EquationTerm::CheckConservativity() {
    bool isConservative = true;
    struct gridcontainer *gc;

    for (len_t i = 0; (gc=GetNextGrid(i)) != nullptr; i++) {
        if (!CheckConservativity(gc->grid)) {
            this->PrintError("%s is not conservative on grid '%s'.", this->name.c_str(), gc->name.c_str());
            isConservative = false;
        }

        delete gc;
    }

    return isConservative;
}

/**
 * Check if the discretization represented by the matrix
 * 'mat' conserves mass.
 *
 * mat: Pre-built matrix representing the discretization to test.
 * rg:  Grid used for the discretization.
 * tol: Relative tolerance to require for agreement.
 */
bool EquationTerm::IsConservative(DREAM::FVM::Matrix *mat, DREAM::FVM::RadialGrid *rg, const real_t tol) {
    Vec sum;
    const len_t n = mat->GetNRows();

    VecCreateSeq(PETSC_COMM_WORLD, n, &sum);

    VecAssemblyBegin(sum);
    VecAssemblyEnd(sum);

    // Compute sum of each row in matrix
    MatGetRowSum(mat->mat(), sum);

    // Fetch values
    PetscScalar y[n];
    PetscInt idx[n];

    for (len_t i = 0; i < n; i++)
        idx[i] = i;

    VecGetValues(sum, n, idx, y);

    // Compute density
    real_t *const* Vp = rg->GetVp();
    real_t I = 0, s = 0;
    len_t offset = 0;
    for (len_t ir = 0; ir < rg->GetNr(); ir++) {
        auto *mg = rg->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                I += y[offset + j*np1 + i] * Vp[ir][j*np1 + i];
                s += Vp[ir][j*np1 + i];
            }
        }

        offset += np1*np2;
    }

    VecDestroy(&sum);

    // We expect I = 0, and s is the total
    // density of a distribution function that
    // is one everywhere (which is the relevant
    // comparison)
    real_t Delta = fabs(I/s);

    if (abs(Delta) >= tol) {
        this->PrintError("I = %e (%f eps)", Delta, Delta/std::numeric_limits<real_t>::epsilon());
        return false;
    } else
        return true;
}
