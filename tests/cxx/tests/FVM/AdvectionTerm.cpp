/**
 * Implementation of tests for the 'AdvectionTerm' class
 * in the DREAM FVM library.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/AdvectionInterpolationCoefficient.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"
#include "FVM/Matrix.hpp"
#include "AdvectionTerm.hpp"
#include "GeneralAdvectionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check if the implementation of the advection term
 * preserves density.
 */
bool AdvectionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;
    GeneralAdvectionTerm *gat = new GeneralAdvectionTerm(grid);
    gat->SetAdvectionBoundaryConditions(
        DREAM::FVM::FLUXGRIDTYPE_P2,
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET,
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET
    );

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = gat->GetNumberOfNonZerosPerRow();
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    for (len_t i = 0; i < 4; i++) {
        // We build the operator in reverse order to avoid causing PETSc
        // to allocate new memory
        const len_t I = 3-i;

        gat->Rebuild(I, 0, nullptr);
        gat->SetMatrixElements(mat, nullptr);
        mat->Assemble();

        const real_t TOLERANCE = 50*NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (I==0?"r" : (I==1?"p1" : (I==2?"p2":"every")));
            this->PrintError("Advection term is not conservative in '%s' component.", dim);

            isConservative = false;
        }
        
        mat->Zero();
    }

    delete mat;
    delete gat;

    return isConservative;
}

/**
 * Check if the implementation of the advection term
 * is correct.
 */
bool AdvectionTerm::CheckValue(DREAM::FVM::Grid *grid) {
    bool isCorrect = true;
    GeneralAdvectionTerm *gat = new GeneralAdvectionTerm(grid);
    gat->SetAdvectionBoundaryConditions(
        DREAM::FVM::FLUXGRIDTYPE_P2,
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET,
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET
    );
    gat->SetAdvectionInterpolationMethod(
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED, 
        DREAM::OptionConstants::AD_INTERP_JACOBIAN_LINEAR, 
        DREAM::FVM::FLUXGRIDTYPE_RADIAL, 0, 1.0);
    gat->SetAdvectionInterpolationMethod(
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED, 
        DREAM::OptionConstants::AD_INTERP_JACOBIAN_LINEAR, 
        DREAM::FVM::FLUXGRIDTYPE_P1, 0, 1.0);
    gat->SetAdvectionInterpolationMethod(
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_INTERP_CENTRED, 
        DREAM::OptionConstants::AD_INTERP_JACOBIAN_LINEAR, 
        DREAM::FVM::FLUXGRIDTYPE_P2, 0, 1.0);

    const len_t ncells = grid->GetNCells();
    real_t *rvec = new real_t[ncells];
    real_t *tvec = new real_t[ncells];
    real_t *f   = new real_t[ncells];

    // XXX here we explicitly assume that the momentum grid
    // is the same att all radii
    const len_t
        nr    = grid->GetNr(),
        n1    = grid->GetMomentumGrid(0)->GetNp1(),
        n2    = grid->GetMomentumGrid(0)->GetNp2();
    const real_t
        *dr   = grid->GetRadialGrid()->GetDr(),
        *d1   = grid->GetMomentumGrid(0)->GetDp1(),
        *d2   = grid->GetMomentumGrid(0)->GetDp2(),
        // Jacobians
        *const* Vp    = grid->GetVp(),
        *const* Vp_fr = grid->GetVp_fr(),
        *const* Vp_f1 = grid->GetVp_f1(),
        *const* Vp_f2 = grid->GetVp_f2();

    // Initialize input and output vectors
    for (len_t i = 0; i < ncells; i++)
        f[i] = 1;

    for (len_t i = 0; i < 4; i++) {
        // Reset vectors
        for (len_t i = 0; i < ncells; i++)
            rvec[i] = tvec[i] = 0;

        gat->Rebuild(i, 0, nullptr);
        gat->SetVectorElements(rvec, f);

        // Evaluate test term
        this->EvaluateFVM(
            nr, n1, n2,
            dr, d1, d2,
            // Jacobians
            [&Vp,   n1](len_t ir, len_t i, len_t j) { return Vp[ir][j*n1+i]; },
            [&Vp_fr,n1](len_t ir, len_t i, len_t j) { return Vp_fr[ir][j*n1+i]; },
            [&Vp_f1,n1](len_t ir, len_t i, len_t j) { return Vp_f1[ir][j*(n1+1)+i]; },
            [&Vp_f2,n1](len_t ir, len_t i, len_t j) { return Vp_f2[ir][j*n1+i]; },
            // Fluxes
            // Phi^(r)
            [&gat,&f,&nr,&n1,&n2](len_t ir, len_t i, len_t j) {
                if (ir == 0 || ir == nr) return 0.0;
                else return 0.5*gat->Fr(ir,i,j)*(f[(ir*n2+j)*n1+i] + f[((ir-1)*n2+j)*n1+i]);
            },
            // Phi^(1)
            [&gat,&f,&n1,&n2](len_t ir, len_t i, len_t j) {
                if (i == 0 || i == n1) return 0.0;
                else return 0.5*gat->F1(ir,i,j)*(f[(ir*n2+j)*n1+i] + f[(ir*n2+j)*n1+(i-1)]);
            },
            // Phi^(2)
            [&gat,&f,&n1,&n2](len_t ir, len_t i, len_t j) {
                if (j == 0 || j == n2) return 0.0;
                else return 0.5*gat->F2(ir,i,j)*(f[(ir*n2+j)*n1+i] + f[(ir*n2+(j-1))*n1+i]);
            },
            // Output vector
            tvec
        );

        // Compare results
        const real_t TOLERANCE = (i+1)*2000*std::numeric_limits<real_t>::epsilon();
        const char coeffnames[4][10] = { "Fr", "F1", "F2", "Fr,F1,F2" };

        for (len_t ir = 0; ir < nr; ir++) {
            for (len_t i2 = 0; i2 < n2; i2++) {
                for (len_t i1 = 0; i1 < n1; i1++) {
                    const len_t idx = (ir*n2 + i2)*n1 + i1;

                    real_t Delta;
                    if (tvec[idx] == 0)
                        Delta = fabs(rvec[idx]);
                    else
                        Delta = fabs(rvec[idx]/tvec[idx]-1);

                    if (Delta > TOLERANCE) {
                        this->PrintError(
                            "'AdvectionTerm' and test vectors are not equal at "
                            "(ir, i2, i1) = "
                            "(" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ") "
                            "with %s =/= 0. Delta = %e.",
                            ir, i2, i1,
                            coeffnames[i],
                            Delta
                        );

                        isCorrect = false;
                        goto ERROR_OCCURED;
                    }
                }
            }
        }
    }

ERROR_OCCURED:

    delete [] tvec;
    delete [] rvec;
    delete [] f;
    delete gat;

    return isCorrect;
}

/**
 * Run all tests for this module.
 */
bool AdvectionTerm::Run(bool) {
    bool success = true;

    // Check the conservativity of the operator on
    // ALL implemented grid combinations
    // (see 'EquationTerm.cpp' for implementation of this
    // routine; see 'UnitTest.cpp' for implementation of
    // the various grid combinations)
    if (!this->EquationTerm::CheckConservativity()) {
        this->PrintError("Conservativity test failed");
        success = false;
    } else
        this->PrintOK("The general advection term conserves density");

    // Check that the term is evaluated correctly
    if (!this->EquationTerm::CheckValue()) {
        this->PrintError("Evaluation test failed");
        success = false;
    } else
        this->PrintOK("The general advection term is evaluated correctly");

    return success;
}

