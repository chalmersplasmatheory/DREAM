/**
 * Test of the 'PXiExternalKineticKinetic' boundary condition which is applied to the
 * boundary between two distribution functions.
 */

#include <iostream>
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
#include "FVM/Equation/BoundaryConditions/PXiExternalLoss.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Matrix.hpp"
#include "GeneralAdvectionTerm.hpp"
#include "GeneralDiffusionTerm.hpp"
#include "PXiExternalKineticKinetic.hpp"


using namespace DREAMTESTS::FVM;
using namespace std;


/**
 * Check that the 'PXiExternalKineticKinetic' boundary condition is internally
 * consistent, i.e. that the particles leaving the hot-tail grid all
 * enter the runaway grid and runaway number.
 */
bool PXiExternalKineticKinetic::CheckConsistency() {
    bool success = true;

    // With same nxi for both RE and hot-tail grids
    if (!Check(&PXiExternalKineticKinetic::CheckConservativity)) {
        this->PrintError("With SAME nxi for both hot-tail and runaway grids.");
        success = false;
    } else
        this->PrintOK("Particle number conserved when nxi is same on both hot-tail and runaway grids.");

    // With different nxi for both RE and hot-tail grids
    if (!Check(&PXiExternalKineticKinetic::CheckConservativity, 34)) {
        this->PrintError("With DIFFERENT nxi for hot-tail and runaway grids.");
        success = false;
    } else
        this->PrintOK("Particle number conserved when nxi is different on hot-tail and runaway grids.");

    // With different nxi AND weird pMax for runaway grid
    if (!Check(&PXiExternalKineticKinetic::CheckConservativity, 42, false)) {
        this->PrintError("With DIFFERENT nxi and pmax for hot-tail and runaway grids.");
        success = false;
    } else
        this->PrintOK("Particle number conserved when nxi and pmax are different on hot-tail and runaway grids.");

    return success;
}

/**
 * Compare the 'PXiExternalKineticKinetic' to the 'PXiExternalLoss' boundary condition.
 * The boundary conditions generally work differently, but can be compared
 * when using the 'f=0' B.C. type for 'PXiExternalLoss' and setting f_RE=0 for
 * 'PXiExternalKineticKinetic'.
 */
bool PXiExternalKineticKinetic::CompareToPXiExternalLoss() {
    return Check(&PXiExternalKineticKinetic::CheckPXiExternalLoss, 0, true);
}

/**
 * General method which sets up the grids and advection/diffusion operators
 * necessary for benchmarking the 'PXiExternalKineticKinetic' boundary condition.
 *
 * checkFunction: Function to call for verifying that the boundary condition
 *                is correctly implemented.
 * nxi_re:        Number of xi points on runaway grid. If 'nxi_re = 0', it is
 *                automatically set to nxi on the hot-tail grid.
 */
bool PXiExternalKineticKinetic::Check(
    bool (PXiExternalKineticKinetic::*checkFunction)(
        DREAM::FVM::Operator*, const string&,
        DREAM::FVM::Grid*, DREAM::FVM::Grid*, DREAM::FVM::Grid*
    ),
    len_t nxi_re, bool sameSizeRE
) {
    bool success = true;
    const len_t nr = 30, np = 4, nxi = 10;

	// This factor is used to scale 'pmax' on the runaway grid. It is useful
	// to make (pmax-pmin)^{RE} different from (pmax-pmin)^{hot} to discover
	// potential implementation errors, but due to how we compare the B.C.
	// to the PXiExternalLoss B.C. we must sometimes make sure that the two
	// grids are of the same size. Otherwise the 'dp' used for diffusion
	// terms on the two grids differ and will cause an error (even if the
	// implementation is actually correct).
	real_t pmaxRE_factor = (sameSizeRE ? 2 : 7);

    if (nxi_re == 0)
        nxi_re = nxi;
    else if (nxi_re == nxi)
        this->PrintWarning("nxi_re = nxi when it probably shouldn't be. nxi_re = " LEN_T_PRINTF_FMT, nxi_re);

    DREAM::FVM::Grid *hottailGrid = this->InitializeGridRCylPXi(nr, np, nxi);
    DREAM::FVM::Grid *runawayGrid = this->InitializeGridRCylPXi(
        nr, np, nxi_re, hottailGrid->GetRadialGrid()->GetBmin(0),
        hottailGrid->GetMomentumGrid(0)->GetP1_f(np),       // pmin
        hottailGrid->GetMomentumGrid(0)->GetP1_f(np)*pmaxRE_factor // pmax
    );
    DREAM::FVM::Grid *fluidGrid   = this->InitializeFluidGrid(nr);

    // Only advection term
    DREAM::FVM::Operator *eqn = new DREAM::FVM::Operator(hottailGrid);
    eqn->AddTerm(new GeneralAdvectionTerm(hottailGrid));

    eqn->RebuildTerms(1, 0, nullptr);   // 1 = F1
    success = success && (this->*checkFunction)(eqn, "Fp", hottailGrid, runawayGrid, fluidGrid);

    delete eqn;
    
    // Only diffusion term
	eqn = new DREAM::FVM::Operator(hottailGrid);
	eqn->AddTerm(new GeneralDiffusionTerm(hottailGrid));

	eqn->RebuildTerms(1, 0, nullptr);	// 1 = Dpp
	success = success && (this->*checkFunction)(eqn, "Dpp", hottailGrid, runawayGrid, fluidGrid);

	delete eqn;

    return success;
}

/**
 * Check that the 'PXiExternalKineticKinetic' B.C. agrees with the 'PXiExternalLoss' B.C.
 * for the given FVM::Operator.
 *
 * eqn:         Operator applied to the hot-tail grid.
 * coeffName:   Name of the coefficient which is set to non-zero.
 * hottailGrid: Grid used for the hot-tail distribution function.
 * runawayGrid: Grid used for the runaway distribution function.
 * fluidGrid:   Fluid grid.
 */
bool PXiExternalKineticKinetic::CheckPXiExternalLoss(
    DREAM::FVM::Operator *eqn, const string& coeffName,
    DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid,
    DREAM::FVM::Grid*
) {
    bool success = true;

    const len_t
        N_hot = hottailGrid->GetNCells(),
        N_re  = runawayGrid->GetNCells();

    DREAM::FVM::BC::PXiExternalLoss *loss =
        new DREAM::FVM::BC::PXiExternalLoss(
            hottailGrid, eqn, 0, 2, nullptr, DREAM::FVM::BC::PXiExternalLoss::BOUNDARY_KINETIC, DREAM::FVM::BC::PXiExternalLoss::BC_F_0
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *cross =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            hottailGrid, hottailGrid, runawayGrid,
            eqn, 0, 1, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
        );

    // Construct test functions
    real_t *f_hot = new real_t[N_hot+N_re];
    real_t *f_re  = f_hot+N_hot;
    for (len_t i = 0; i < N_hot; i++)
        f_hot[i] = 1.0 + i;
    for (len_t i = 0; i < N_re; i++)
        f_re[i] = 0;

    DREAM::FVM::Matrix *lossMat = new DREAM::FVM::Matrix(N_hot, N_hot, 5);
    DREAM::FVM::BlockMatrix *crossMat = new DREAM::FVM::BlockMatrix();

    crossMat->CreateSubEquation(N_hot, 5, 0);
    crossMat->CreateSubEquation(N_re, 15, 1);
    crossMat->ConstructSystem();

    loss->AddToMatrixElements(lossMat, nullptr);
    cross->AddToMatrixElements(crossMat, nullptr);

    lossMat->Assemble();
    crossMat->Assemble();

    real_t *PhiLoss  = lossMat->Multiply(N_hot, f_hot);
    real_t *PhiCross = crossMat->Multiply(N_hot+N_re, f_hot);

    // Compare fluxes
    const real_t TOLERANCE = 100 * std::numeric_limits<real_t>::epsilon();
    for (len_t i = 0; i < N_hot; i++) {
        real_t Delta;
        if (PhiLoss[i] == 0)
            Delta = abs(PhiCross[i]);
        else
            Delta = abs(PhiCross[i] / PhiLoss[i] - 1.0);

        if (Delta >= TOLERANCE) {
            this->PrintError(
                "PXiExternalKineticKinetic does not agree with PXiExternalLoss at index "
                LEN_T_PRINTF_FMT " with %s =/= 0. Delta = %e.",
                i, coeffName.c_str(), Delta
            );

            success = false;
            break;
        }
    }

    delete [] f_hot;
    
    delete [] PhiLoss;
    delete [] PhiCross;

    delete lossMat;
    delete crossMat;
    delete loss;
    delete cross;

    return success;
}

/**
 * Check that the number of particles entering the fluid/runaway grids
 * is the same as the number of particles leaving the hot-tail grid.
 */
bool PXiExternalKineticKinetic::CheckConservativity(
    DREAM::FVM::Operator *eqn, const string& coeffName,
    DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid,
    DREAM::FVM::Grid *fluidGrid
) {
    const real_t TOLERANCE = 100 * std::numeric_limits<real_t>::epsilon();
    bool success = true;

    const len_t
        N_hot  = hottailGrid->GetNCells(),
        N_re   = runawayGrid->GetNCells(),
        N_dens = fluidGrid->GetNCells(),
        N_hot_mom = hottailGrid->GetMomentumGrid(0)->GetNCells(),
        N_re_mom  = runawayGrid->GetMomentumGrid(0)->GetNCells();

    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();
    uqh->InsertUnknown("f_hot", "0", hottailGrid);
    uqh->InsertUnknown("f_re", "0", runawayGrid);
    uqh->InsertUnknown("n_re", "0", fluidGrid);

    const len_t id_f_hot = uqh->GetUnknownID("f_hot");
    const len_t id_f_re  = uqh->GetUnknownID("f_re");
    const len_t id_n_re  = uqh->GetUnknownID("n_re");

    // Create boundary conditions
    DREAM::FVM::BC::PXiExternalKineticKinetic *lower =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            hottailGrid, hottailGrid, runawayGrid,
            eqn, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *upper =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            runawayGrid, hottailGrid, runawayGrid,
            eqn, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *density =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            fluidGrid, hottailGrid, runawayGrid,
            eqn, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_DENSITY
        );

    // Construct test functions
    real_t *f_hot = new real_t[N_hot+N_re+N_dens];
    real_t *f_re  = f_hot+N_hot;
    real_t *n_re  = f_re+N_re;
    for (len_t i = 0; i < N_hot; i++)
        f_hot[i] = 1.0 + i;
		//f_hot[i] = 1.0;
    for (len_t i = 0; i < N_re; i++)
        //f_re[i]  = f_hot[N_hot-1] - i*(real_t(N_hot)/real_t(N_re));
        f_re[i] = 1.0 + i*i;
		//f_re[i] = 1.0;
    for (len_t i = 0; i < N_dens; i++)
        n_re[i] = 0;

    DREAM::FVM::BlockMatrix *lowerMat = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *upperMat = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *densMat  = new DREAM::FVM::BlockMatrix();

    lowerMat->CreateSubEquation(N_hot, 5, id_f_hot);
    lowerMat->CreateSubEquation(N_re, 5, id_f_re);
    lowerMat->ConstructSystem();

    upperMat->CreateSubEquation(N_hot, 5, id_f_hot);
    upperMat->CreateSubEquation(N_re, 5, id_f_re);
    upperMat->ConstructSystem();

    densMat->CreateSubEquation(N_hot, 1, id_f_hot);
    densMat->CreateSubEquation(N_re, 1, id_f_re);
    densMat->CreateSubEquation(N_dens, N_hot_mom+N_re_mom, id_n_re);
    densMat->ConstructSystem();

    lowerMat->SetOffset(0, 0);
    upperMat->SetOffset(N_hot, N_hot);
    densMat->SetOffset(N_hot+N_re, N_hot+N_re);

    lower->AddToMatrixElements(lowerMat, nullptr);
    upper->AddToMatrixElements(upperMat, nullptr);
    density->AddToMatrixElements(densMat, nullptr);

    lowerMat->Assemble();
    upperMat->Assemble();
    densMat->Assemble();

    real_t *PhiLower  = lowerMat->Multiply(N_hot+N_re, f_hot);
    real_t *PhiUpper_ = upperMat->Multiply(N_hot+N_re, f_hot);
    real_t *PhiDens_  = densMat->Multiply(N_hot+N_re+N_dens, f_hot);

    real_t *PhiUpper = PhiUpper_ + N_hot;
    real_t *PhiDens  = PhiDens_  + N_hot+N_re;

    lowerMat->View(DREAM::FVM::Matrix::BINARY_MATLAB, "petsc_lower");
    upperMat->View(DREAM::FVM::Matrix::BINARY_MATLAB, "petsc_upper");
    densMat->View(DREAM::FVM::Matrix::BINARY_MATLAB, "petsc_dens");

    // Compare lower/upper fluxes element-by-element
    //   PhiL * VpL * dxiL = sum[ PhiU * VpU * dxiBar ]
    real_t *PhiUpper_conv = ConvertFlux(PhiUpper, runawayGrid, hottailGrid);
	if ((N_hot < N_re && N_hot%N_re == 0) || (N_hot > N_re && N_re%N_hot == 0)) {
		for (len_t i = 0; i < N_hot; i++) {
			real_t Delta;

			if (PhiLower[i] == 0)
				Delta = abs(PhiUpper_conv[i]);
			else
				Delta = abs(PhiUpper_conv[i] / PhiLower[i] + 1.0);

			if (Delta > TOLERANCE) {
				// XXX here we assume that all momentum grids are the same
				len_t
					ir  = i / N_hot_mom,
					ixi = (i-ir*N_hot_mom) / (hottailGrid->GetMomentumGrid(0)->GetNp1());

				this->PrintError(
					"hot -> re: deviation in element comparison at i = " LEN_T_PRINTF_FMT " "
					"(ir, ixi = " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ") "
					"with %s =/= 0. Delta = %e.",
					i, ir, ixi, coeffName.c_str(), Delta
				);

				success = false;
				break;
			}
		}
	}

    // Integrate Phi_hot and Phi_RE over momentum (p and xi)
    real_t *lowerI = hottailGrid->IntegralMomentum(PhiLower);
    real_t *upperI = runawayGrid->IntegralMomentum(PhiUpper);
	real_t *upperI_conv = hottailGrid->IntegralMomentum(PhiUpper_conv);

    // Compare integrated fluxes
    for (len_t ir = 0; ir < N_dens; ir++) {
        real_t Delta;

        // f_hot  -->  n_re
        if (lowerI[ir] == 0)
            Delta = abs(PhiDens[ir]);
        else
            Delta = abs(PhiDens[ir] / lowerI[ir] + 1);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "hot -> fluid: deviation at ir = " LEN_T_PRINTF_FMT " "
                "with %s =/= 0. Delta = %e.",
                ir, coeffName.c_str(), Delta
            );

            success = false;
            break;
        }

        // f_hot  -->  f_re
        if (lowerI[ir] == 0)
            Delta = abs(upperI[ir]);
        else
            Delta = abs(upperI[ir] / lowerI[ir] + 1);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "hot -> re: deviation at ir = " LEN_T_PRINTF_FMT " "
                "with %s =/= 0. Delta = %e.",
                ir, coeffName.c_str(), Delta
            );
            
            success = false;
            break;
        }

		// f_hot  -->  f_re
		// (calculated by converting Phi^{RE} to Phi^{hot})
		if (lowerI[ir] == 0)
			Delta = abs(upperI_conv[ir]);
		else
			Delta = abs(upperI_conv[ir] / lowerI[ir] + 1);

		if (Delta > TOLERANCE) {
			this->PrintError(
				"hot -> re [conv]: deviation at ir = " LEN_T_PRINTF_FMT " "
				"with %s =/= 0. Delta = %e.",
				ir, coeffName.c_str(), Delta
			);
		}
    }

    delete [] lowerI;
    delete [] upperI;

    delete [] PhiUpper_conv;
    delete [] PhiDens_;
    delete [] PhiUpper_;
    delete [] PhiLower;

    delete [] f_hot;

    delete density;
    delete upper;
    delete lower;

    delete uqh;

    return success;
}

/**
 * Convert a flux evaluated on 'grid1' to a flux on 'grid2'. This
 * can be used to compare fluxes on different, but adjacent, grids.
 */
real_t *PXiExternalKineticKinetic::ConvertFlux(
    const real_t *Phi1, DREAM::FVM::Grid *grid1,
    DREAM::FVM::Grid *grid2
) {
    len_t offset1 = 0, offset2 = 0;
    const len_t
        nr   = grid2->GetNr(),
        N2   = grid2->GetNCells();

    real_t *Phi2 = new real_t[N2];
    for (len_t i = 0; i < N2; i++)
        Phi2[i] = 0;

    for (len_t ir = 0; ir < nr; ir++) {
        DREAM::FVM::MomentumGrid *mg1 = grid1->GetMomentumGrid(ir);
        DREAM::FVM::MomentumGrid *mg2 = grid2->GetMomentumGrid(ir);

        const len_t
            np1  = mg1->GetNp1(),
            nxi1 = mg1->GetNp2(),
            np2  = mg2->GetNp1(),
            nxi2 = mg2->GetNp2();

        const real_t
            *Vp1   = grid1->GetVp(ir),
            *Vp2   = grid2->GetVp(ir),
            *dp1   = mg1->GetDp1(),
            *dp2   = mg2->GetDp1(),
            *dxi2  = mg2->GetDp2(),
            *xi1_f = mg1->GetP2_f(),
            *xi2_f = mg2->GetP2_f();

        for (len_t j = 0; j < nxi2; j++) {
            len_t idx2   = j*np2 + np2-1;

            //////////////////
            // SUM OVER J
            //////////////////
            #define OVERLAPPING(j,J) (xi1_f[(J)+1]>xi2_f[(j)]) && (xi1_f[(J)]<xi2_f[(j)+1])
            // Locate starting xi
            len_t J = 0;
            while (J < nxi1 && not OVERLAPPING(j,J))
                J++;
            
            real_t s = 0;
            while (J < nxi1 && OVERLAPPING(j,J)) {
                len_t idx1 = J*np1;
                real_t dxiBar = min(xi1_f[J+1], xi2_f[j+1]) - max(xi1_f[J], xi2_f[j]);

                s += Phi1[offset1+idx1] * Vp1[idx1] * dxiBar * dp1[0];

                J++;
            }
            #undef OVERLAPPING

            Phi2[offset2+idx2] = s / (Vp2[idx2]*dxi2[j]*dp2[np2-1]);
        }

        offset1 += np1*nxi1;
        offset2 += np2*nxi2;
    }

    return Phi2;
}


/**
 * Run this test.
 */
bool PXiExternalKineticKinetic::Run(bool) {
    bool success = true;

    // Compare to the PXiExternalLoss boundary condition
    if (!CompareToPXiExternalLoss()) {
        this->PrintError("Flux does not agree with the PXiExternalLoss boundary condition.");
        success = false;
    } else
        this->PrintOK("Flux agrees with PXiExternalLoss boundary condition.");

    // Check that B.C. is conservative on hot-tail, runaway and fluid grids.
    // Runaway grid should be tested with different number of xi points.
    if (!CheckConsistency())
        success = false;
    else
        this->PrintOK("Boundary condition is internally consistent.");

    return success;
}

