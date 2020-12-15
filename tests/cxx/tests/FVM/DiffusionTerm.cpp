/**
 * Implementation of tests for the 'DiffusionTerm' class
 * in the DREAM FVM library.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Equation/BoundaryConditions/PXiInternalTrapping.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Matrix.hpp"
#include "DiffusionTerm.hpp"
#include "GeneralDiffusionTerm.hpp"


using namespace DREAMTESTS::FVM;

/**
 * Check if the implementation of the diffusion term
 * preserves density.
 */
bool DiffusionTerm::CheckConservativity(DREAM::FVM::Grid *grid) {
    bool isConservative = true;
    DREAM::FVM::Operator *Op = new DREAM::FVM::Operator(grid);
    Op->AddTerm(new GeneralDiffusionTerm(grid,1.0));
    if(grid->HasTrapped())
        Op->AddBoundaryCondition(new DREAM::FVM::BC::PXiInternalTrapping(grid,Op));

    const len_t ncells = grid->GetNCells();
    const len_t NNZ_PER_ROW = Op->GetNumberOfNonZerosPerRow();
    DREAM::FVM::Matrix *mat = new DREAM::FVM::Matrix(ncells, ncells, NNZ_PER_ROW);

    for (len_t i = 0; i < 6; i++) {
        Op->RebuildTerms(5-i, 0, nullptr);
        Op->SetMatrixElements(mat, nullptr);
        mat->Assemble();

        const real_t TOLERANCE = 10*NNZ_PER_ROW*ncells * std::numeric_limits<real_t>::epsilon();

        if (!IsConservative(mat, grid, TOLERANCE)) {
            const char *dim = (i==0?"rr" : (i==1?"11" : (i==2?"22": (i==3?"12" : (i==4?"21":"every")))));
            this->PrintError("Diffusion term is not conservative in '%s' component.", dim);

            isConservative = false;
        }
        if((5-i)==2)
            mat->View(DREAM::FVM::Matrix::BINARY_MATLAB);
        mat->Zero();
    }

    delete mat;
    delete Op;

    return isConservative;
}

/**
 * Check if the implementation of the advection term
 * is correct.
 */
bool DiffusionTerm::CheckValue(DREAM::FVM::Grid *grid) {
    bool isCorrect = true;
    GeneralDiffusionTerm *gdt = new GeneralDiffusionTerm(grid);

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
        *dr_f = grid->GetRadialGrid()->GetDr_f(),
        *d1   = grid->GetMomentumGrid(0)->GetDp1(),
        *d1_f = grid->GetMomentumGrid(0)->GetDp1_f(),
        *d2   = grid->GetMomentumGrid(0)->GetDp2(),
        *d2_f = grid->GetMomentumGrid(0)->GetDp2_f(),
        // Jacobians
        *const* Vp    = grid->GetVp(),
        *const* Vp_fr = grid->GetVp_fr(),
        *const* Vp_f1 = grid->GetVp_f1(),
        *const* Vp_f2 = grid->GetVp_f2();

    // Initialize input and output vectors
    for (len_t i = 0; i < ncells; i++)
        f[i] = 1.0 + i;

    for (len_t i = 0; i < 6; i++) {
        // Reset vectors
        for (len_t i = 0; i < ncells; i++)
            rvec[i] = tvec[i] = 0;

        gdt->Rebuild(i, 0, nullptr);
        gdt->SetVectorElements(rvec, f);

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
            [&gdt,&f,&dr_f,&nr,&n1,&n2](len_t ir, len_t i, len_t j) {
                if (ir == 0 || ir == nr) return 0.0;
                else return -gdt->Drr(ir,i,j) * (
                    f[(ir*n2+j)*n1+i] - f[((ir-1)*n2+j)*n1+i]
                ) / dr_f[ir-1];
            },
            // Phi^(1)
            [&gdt,&f,&d1_f,&d2_f,&n1,&n2](len_t ir, len_t i, len_t j) {
                real_t S = 0;

                if (i > 0 && i < n1) {
                    S -= gdt->D11(ir,i,j) * (
                        f[(ir*n2+j)*n1+i] - f[(ir*n2+j)*n1+(i-1)]
                    ) / d1_f[i-1];

                    if (j > 0 && j < n2-1)
                        S -= gdt->D12(ir,i,j) * (
                            f[(ir*n2+(j+1))*n1 + i] + f[(ir*n2+(j+1))*n1 + (i-1)] -
                            f[(ir*n2+(j-1))*n1 + i] - f[(ir*n2+(j-1))*n1 + (i-1)]
                        ) / (d2_f[j]+d2_f[j-1]);
                }

                return S;
            },
            // Phi^(2)
            [&gdt,&f,&d1_f,&d2_f,&n1,&n2](len_t ir, len_t i, len_t j) {
                real_t S = 0;

                if (j > 0 && j < n2) {
                    S -= gdt->D22(ir,i,j) * (
                        f[(ir*n2+j)*n1+i] - f[(ir*n2+(j-1))*n1+i]
                    ) / d2_f[j-1];

                    if (i > 0 && i < n1-1)
                        S -= gdt->D21(ir,i,j) * (
                            f[(ir*n2+j)*n1 + (i+1)] + f[(ir*n2+(j-1))*n1 + (i+1)] -
                            f[(ir*n2+j)*n1 + (i-1)] - f[(ir*n2+(j-1))*n1 + (i-1)]
                        ) / (d1_f[i] + d1_f[i-1]);
                }
                
                return S;
            },
            // Output vector
            tvec
        );

        // Compare results
        const real_t TOLERANCE = sqrt(std::numeric_limits<real_t>::epsilon());
        const char coeffnames[6][18] = { "Drr", "D11", "D22", "D12", "D21", "all coefficients" };

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
                            "'DiffusionTerm' and test vectors are not equal at "
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
    delete gdt;

    return isCorrect;
}

/**
 * Run all tests for this module.
 */
bool DiffusionTerm::Run(bool) {
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
        this->PrintOK("The general diffusion term conserves density");

    // Check that the term is evaluated correctly
    if (!this->EquationTerm::CheckValue()) {
        this->PrintError("Evaluation test failed");
        success = false;
    } else
        this->PrintOK("The general diffusion term is evaluated correctly");

    return success;
}

