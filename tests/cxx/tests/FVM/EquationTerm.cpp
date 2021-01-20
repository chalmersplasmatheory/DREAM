/**
 * Implementation of general tests for all 'EquationTerm' classes.
 */

#include <functional>
#include <petsc.h>
#include "FVM/Matrix.hpp"
#include "FVM/Grid/Grid.hpp"
#include "EquationTerm.hpp"


using namespace DREAMTESTS::FVM;
using namespace std;


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
 * Check that this equation term is evaluated correctly.
 */
bool EquationTerm::CheckValue() {
    bool isCorrect = true;
    struct gridcontainer *gc;

    for (len_t i = 0; (gc=GetNextGrid(i)) != nullptr; i++) {
        if(gc->grid->HasTrapped())
            continue; // CHECK NOT IMPLEMENTED FOR INHOMOGENEOUS FIELDS 
        if (!CheckValue(gc->grid)) {
            this->PrintError("%s is not evaluated correctly on grid '%s'.", this->name.c_str(), gc->name.c_str());
            isCorrect = false;
        }

        delete gc;
    }

    return isCorrect;
}

/**
 * Evaluate the finite volume method with the given fluxes
 * on the specified grid.
 *
 * nr, n1, n2:          Grid sizes.
 * Vp:                  Grid jacobian.
 * Vp_fr, Vp_f1, Vp_f2: Jacobians evaluated on flux grids.
 * fluxR, flux1, flux2: Phase space flux functions.
 * vec:                 Vector containing evaluated equation term
 *                      on return.
 */
void EquationTerm::EvaluateFVM(
    const len_t nr, const len_t n1, const len_t n2,
    const real_t *dr, const real_t *d1, const real_t *d2,
    function<real_t(len_t, len_t, len_t)> Vp,
    function<real_t(len_t, len_t, len_t)> Vp_fr,
    function<real_t(len_t, len_t, len_t)> Vp_f1,
    function<real_t(len_t, len_t, len_t)> Vp_f2,
    function<real_t(len_t, len_t, len_t)> fluxR,
    function<real_t(len_t, len_t, len_t)> flux1,
    function<real_t(len_t, len_t, len_t)> flux2,
    real_t *vec
) {
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                const len_t idx = (ir*n2 + j)*n1 + i;

                vec[idx] =
                    // r-term
                    (Vp_fr(ir+1,i,j)*fluxR(ir+1,i,j) - Vp_fr(ir,i,j)*fluxR(ir,i,j))
                    / (Vp(ir,i,j)*dr[ir]) +

                    // p1-term
                    (Vp_f1(ir,i+1,j)*flux1(ir,i+1,j) - Vp_f1(ir,i,j)*flux1(ir,i,j))
                    / (Vp(ir,i,j)*d1[i]) +

                    // p2-term
                    (Vp_f2(ir,i,j+1)*flux2(ir,i,j+1) - Vp_f2(ir,i,j)*flux2(ir,i,j))
                    / (Vp(ir,i,j)*d2[j]);
            }
        }
    }
}


/**
 * Check if the discretization represented by the matrix
 * 'mat' conserves mass.
 *
 * mat:  Pre-built matrix representing the discretization to test.
 * grid: Grid used for the discretization.
 * tol:  Relative tolerance to require for agreement.
 */
bool EquationTerm::IsConservative(DREAM::FVM::Matrix *mat, DREAM::FVM::Grid *grid, const real_t tol) {
    Vec integratedTermVec;
    Vec densityIntegralVec;
    const len_t n = mat->GetNRows();
    const len_t m = mat->GetNCols();
    
    VecCreateSeq(PETSC_COMM_WORLD, m, &densityIntegralVec);


    VecCreateSeq(PETSC_COMM_WORLD, n, &integratedTermVec);
    VecAssemblyBegin(integratedTermVec);
    VecAssemblyEnd(integratedTermVec);

    // Set elements in the density-integration integrannd
    real_t *const* Vp = grid->GetVp();
    len_t offset = 0;
    real_t Volume = 0;
    for (len_t ir = 0; ir < grid->GetNr(); ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        const real_t dr = grid->GetRadialGrid()->GetDr(ir);
        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1; i++) {
                const real_t dp1 = mg->GetDp1(i);
                const real_t dp2 = mg->GetDp2(j);                
                PetscScalar dV = Vp[ir][j*np1+i]*dr*dp1*dp2;
                PetscInt idx = offset + j*np1 + i;
                Volume += dV;
                VecSetValue(densityIntegralVec, idx, dV, INSERT_VALUES);
            }
        offset += np1*np2;
    }
    VecAssemblyBegin(densityIntegralVec);
    VecAssemblyEnd(densityIntegralVec);

    // integratedTerm = (equation term matrix)^T * densityIntegral
    MatMultTranspose(mat->mat(), densityIntegralVec, integratedTermVec);

    // Fetch values
    PetscScalar *integratedTerm = new PetscScalar[n];
    PetscInt *idx = new PetscInt[n];
    for (len_t i = 0; i < n; i++)
        idx[i] = i;

    VecGetValues(integratedTermVec, n, idx, integratedTerm);

    VecDestroy(&densityIntegralVec);
    VecDestroy(&integratedTermVec);
    delete [] idx;
    

    bool success = true;
    offset = 0;
    for (len_t ir = 0; ir < grid->GetNr(); ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1; i++) {
                real_t Delta = fabs(integratedTerm[offset + np1*j + i]) / Volume;
                if(fabs(Delta)>=tol){
                    success = false;
                    if(i==0)
                        this->PrintError(
                            "I = %e (%f eps), " 
                            "in element ir = %u, j = %u", 
                            Delta, Delta/std::numeric_limits<real_t>::epsilon(),
                            ir, j
                        );
                }
            }
        offset += np1*np2;
    }

    delete [] integratedTerm;
    return success;
}

