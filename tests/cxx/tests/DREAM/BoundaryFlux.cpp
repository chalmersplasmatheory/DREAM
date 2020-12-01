/**
 * Tests that the flux across grid boundaries are evaluated
 * consistently.
 */

#include <iostream>
#include "DREAM/Equations/Fluid/DensityFromBoundaryFluxPXI.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "BoundaryFlux.hpp"
#include "../FVM/GeneralAdvectionTerm.hpp"
#include "../FVM/GeneralDiffusionTerm.hpp"


using namespace DREAMTESTS::_DREAM;
using namespace std;


/**
 * Check that the flux out of a P/XI grid, into a fluid
 * grid, is consistent.
 */
bool BoundaryFlux::CheckPXiToFluid() {
    bool success = true;
    const len_t nr = 30, np = 40, nxi = 10;

    DREAM::FVM::Grid *kineticGrid = this->InitializeGridRCylPXi(nr, np, nxi);
    DREAM::FVM::Grid *fluidGrid   = this->InitializeFluidGrid(nr);

    // Only advection term
    const char advCoeffs[4][10] = {"Fr", "F1", "F2", "Fr,F1,F2"};
    DREAM::FVM::Operator *eqn = new DREAM::FVM::Operator(kineticGrid);
    eqn->AddTerm(new FVM::GeneralAdvectionTerm(kineticGrid));

    for (len_t i = 0; i < 4; i++) {
        eqn->RebuildTerms(i, 0, nullptr);

        success = success && CheckPXiToFluid(eqn, advCoeffs[i], kineticGrid, fluidGrid);
    }

    delete eqn;

    // Only diffusion term
    const char diffCoeffs[6][20] = {"Drr", "D11", "D22", "D12", "D21", "all diff. coeffs."};
    eqn = new DREAM::FVM::Operator(kineticGrid);
    eqn->AddTerm(new FVM::GeneralDiffusionTerm(kineticGrid));

    for (len_t i = 0; i < 6; i++) {
        eqn->RebuildTerms(i, 0, nullptr);

        success = success && CheckPXiToFluid(eqn, diffCoeffs[i], kineticGrid, fluidGrid);
    }

    // XXX For some reason this test fails with suspiciously large errors
    // despite all the other tests passing with (relative) errors less than ~1e-12
    //
    // Combined advection & diffusion
    eqn->AddTerm(new FVM::GeneralAdvectionTerm(kineticGrid));
    // 100 = guaranteed set all coefficients non-zero
    eqn->RebuildTerms(100, 0, nullptr);
    success = success && CheckPXiToFluid(eqn, "all coefficients", kineticGrid, fluidGrid);

    delete eqn;

    return success;
}

/**
 * Check that the flux out of a P/XI grid, into a fluid grid,
 * is consistent. Do this for the specified equation.
 *
 */
bool BoundaryFlux::CheckPXiToFluid(
    const DREAM::FVM::Operator *eqn, const string& coeffnames,
    DREAM::FVM::Grid *kineticGrid, DREAM::FVM::Grid *fluidGrid
) {
    bool success = true;
    const len_t
        M = fluidGrid->GetNCells(),
        N = kineticGrid->GetNCells();

    DREAM::DensityFromBoundaryFluxPXI *dens =
        new DREAM::DensityFromBoundaryFluxPXI(fluidGrid, kineticGrid, eqn, 0,0);

    DREAM::FVM::BC::PXiExternalLoss *loss =
        new DREAM::FVM::BC::PXiExternalLoss(kineticGrid, eqn, 0);

    // Construct test function
    real_t *f = new real_t[N];
    for (len_t i = 0; i < N; i++)
        f[i] = 1.0 + i;

    ///////////////////////////////////
    // MATRIX TEST
    ///////////////////////////////////
    // XXX here we explicitly assume that all momentum grids are the same
    DREAM::FVM::Matrix *densMat = new DREAM::FVM::Matrix(M, N, dens->GetNumberOfNonZerosPerRow());
    DREAM::FVM::Matrix *lossMat =
        new DREAM::FVM::Matrix(N, N, kineticGrid->GetMomentumGrid(0)->GetNp2());

    dens->Rebuild(0, 0, nullptr);
    loss->Rebuild(0, nullptr);

    // Construct matrices
    dens->SetMatrixElements(densMat, nullptr);
    loss->AddToMatrixElements(lossMat, nullptr);

    densMat->Assemble();
    lossMat->Assemble();

    // Apply operators to test function
    real_t *densI    = densMat->Multiply(N, f);
    real_t *lossGrid = lossMat->Multiply(N, f);

    // Integrate lost particles
    real_t *lossI = kineticGrid->IntegralMomentum(lossGrid);

    // Compare results for 'SetMatrixElements()' and 'AddToMatrixElements()'
    // in 'DensityFromBoundaryFluxPXI' and 'PXiExternalLoss' respectively
    success = success &&
        VerifyVectorsAreOpposite("Matrix calculations", coeffnames, M, N, densI, lossI);

    ///////////////////////////////////
    // VECTOR TEST
    ///////////////////////////////////
    for (len_t i = 0; i < M; i++)
        densI[i] = 0;
    for (len_t i = 0; i < N; i++)
        lossGrid[i] = 0;

    dens->SetVectorElements(densI, f);
    loss->AddToVectorElements(lossGrid, f);

    // Compare result of 'SetVectorElements()' in 'DensityFromBoundaryFluxPXI'
    // to the previously evaluated 'AddToMatrixElements()' in 'PXiExternalLoss'
    success = success &&
        VerifyVectorsAreOpposite("Matrix/vector calculations", coeffnames, M, N, densI, lossI);

    delete [] lossI;

    lossI = kineticGrid->IntegralMomentum(lossGrid);
    success = success && VerifyVectorsAreOpposite("Vector calculations", coeffnames, M, N, densI, lossI);

    delete [] densI;
    delete [] lossGrid;

    delete lossMat;
    delete densMat;
    delete loss;
    delete dens;

    return success;
}

/**
 * Check that the vectors are equal and opposite.
 */
bool BoundaryFlux::VerifyVectorsAreOpposite(
    const string& calcname, const string& coeffnames,
    const len_t NR, const len_t nCells, 
    const real_t *vec1, const real_t *vec2
) {
    const real_t TOLERANCE = nCells*std::numeric_limits<real_t>::epsilon();

    real_t maxDelta = 0;
    for (len_t ir = 0; ir < NR; ir++) {
        real_t Delta;
        if (vec2[ir] == 0)
            Delta = fabs(vec1[ir]);
        else
            // The fluxes are supposed to cancel (i.e. be equal and opposite)
            Delta = fabs(vec1[ir]/vec2[ir] + 1);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "p/xi -> fluid: %s deviate at ir = " LEN_T_PRINTF_FMT " "
                "with %s =/= 0. Delta = %e.",
                calcname.c_str(), ir, coeffnames.c_str(), Delta
            );

            return false;
        }

        if (Delta > maxDelta)
            maxDelta = Delta;
    }
    
    return true;
}

/**
 * Run this test.
 */
bool BoundaryFlux::Run(bool) {
    bool success = true;

    if (!CheckPXiToFluid()) {
        this->PrintError("p/xi -> fluid: Flux is not evaluated consistently");
        success = false;
    } else
        this->PrintOK("p/xi -> fluid: Flux is evaluated consistently");

    return success;
}

