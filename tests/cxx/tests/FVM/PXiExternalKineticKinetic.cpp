/**
 * Test of the 'PXiExternalKineticKinetic' boundary condition which is applied to the
 * boundary between two distribution functions.
 */

#include <iostream>
#include "FVM/Equation/BoundaryConditions/PXiExternalKineticKinetic.hpp"
//#include "FVM/Equation/BoundaryConditions/PXiExternalKineticLower.hpp"
//#include "FVM/Equation/BoundaryConditions/PXiExternalKineticUpper.hpp"
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

    // With different nxi for both RE and hot-tail grids (nxi_RE > nxi_hot)
    if (!Check(&PXiExternalKineticKinetic::CheckConservativity, 20)) {
        this->PrintError("With MORE nxi for runaway than the hot-tail grid.");
        success = false;
    } else
        this->PrintOK("Particle number conserved when nxi is MORE on runaway than the hot-tail grid.");

    // With nxi_RE < nxi_hot
    if (!Check(&PXiExternalKineticKinetic::CheckConservativity, 5)) {
        this->PrintError("With LESS nxi for runaway than the hot-tail grid.");
        success = false;
    } else
        this->PrintOK("Particle number conserved when nxi is LESS on runaway than hot-tail grid.");

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
 * Compare the 'PXiExternalKineticKinetic' to the reference implementations
 * for the boundary conditions.
 */
/*bool PXiExternalKineticKinetic::CompareToReference() {
    bool success = true;

    // Same nxi
    if (!Check(&PXiExternalKineticKinetic::CheckWithReference, 0)) {
        this->PrintError("With SAME nxi for both hot-tail and runaway grids.");
        success = false;
    }

    // Different nxi
    if (!Check(&PXiExternalKineticKinetic::CheckWithReference, 33)) {
        this->PrintError("With DIFFERENT nxi for both hot-tail and runaway grids.");
        success = false;
    }

    // Different nxi and differen pmax
    if (!Check(&PXiExternalKineticKinetic::CheckWithReference, 39, false)) {
        this->PrintError("With DIFFERENT nxi and pmax for both hot-tail and runaway grids.");
        success = false;
    }

    return success;
}*/

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
	/*eqn = new DREAM::FVM::Operator(hottailGrid);
	eqn->AddTerm(new GeneralDiffusionTerm(hottailGrid));

	eqn->RebuildTerms(1, 0, nullptr);	// 1 = Dpp
	success = success && (this->*checkFunction)(eqn, "Dpp", hottailGrid, runawayGrid, fluidGrid);

	delete eqn;*/

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
 * agrees with the more specific reference implementations.
 */
/*bool PXiExternalKineticKinetic::CheckWithReference(
    DREAM::FVM::Operator *eqn, const string& coeffName,
    DREAM::FVM::Grid *hottailGrid, DREAM::FVM::Grid *runawayGrid,
    DREAM::FVM::Grid* fluidGrid
) {
    bool success = true;

    const len_t
        N_hot = hottailGrid->GetNCells(),
        N_re  = runawayGrid->GetNCells();

    DREAM::FVM::BC::PXiExternalKineticLower *refLow =
        new DREAM::FVM::BC::PXiExternalKineticLower(
            hottailGrid, runawayGrid, eqn, 0, 1
        );
    DREAM::FVM::BC::PXiExternalKineticUpper *refUpp =
        new DREAM::FVM::BC::PXiExternalKineticUpper(
            hottailGrid, runawayGrid, eqn, 0, 1
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *crossLower =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            hottailGrid, hottailGrid, runawayGrid,
            eqn, 0, 1, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *crossUpper =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            runawayGrid, hottailGrid, runawayGrid,
            eqn, 0, 1, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
        );

    // Construct test functions
    real_t *f_hot = new real_t[N_hot+N_re];
    real_t *f_re  = f_hot+N_hot;
    for (len_t i = 0; i < N_hot; i++)
        f_hot[i] = 1.0 + i;
    for (len_t i = 0; i < N_re; i++)
        f_re[i] = N_hot + i;

    // Evaluate operators
    DREAM::FVM::BlockMatrix *refLowMat = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *refUppMat = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *crossLowMat = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *crossUppMat = new DREAM::FVM::BlockMatrix();

    refLowMat->CreateSubEquation(N_hot, 5, 0);
    refLowMat->CreateSubEquation(N_re, 15, 1);
    refLowMat->ConstructSystem();

    refUppMat->CreateSubEquation(N_hot, 5, 0);
    refUppMat->CreateSubEquation(N_re, 15, 1);
    refUppMat->ConstructSystem();

    crossLowMat->CreateSubEquation(N_hot, 5, 0);
    crossLowMat->CreateSubEquation(N_re, 15, 1);
    crossLowMat->ConstructSystem();

    crossUppMat->CreateSubEquation(N_hot, 5, 0);
    crossUppMat->CreateSubEquation(N_re, 15, 1);
    crossUppMat->ConstructSystem();

    refLow->AddToMatrixElements(refLowMat, nullptr);
    refUpp->AddToMatrixElements(refUppMat, nullptr);
    crossLower->AddToMatrixElements(crossLowMat, nullptr);
    crossUpper->AddToMatrixElements(crossUppMat, nullptr);

    refLowMat->Assemble();
    refUppMat->Assemble();
    crossLowMat->Assemble();
    crossUppMat->Assemble();

    real_t *PhiRefLow = refLowMat->Multiply(N_hot+N_re, f_hot);
    real_t *PhiRefUpp = refUppMat->Multiply(N_hot+N_re, f_hot);
    real_t *PhiKinLow = crossLowMat->Multiply(N_hot+N_re, f_hot);
    real_t *PhiKinUpp = crossUppMat->Multiply(N_hot+N_re, f_hot);

    // Compare fluxes (lower grid)
    const real_t TOLERANCE = 100*std::numeric_limits<real_t>::epsilon();
    for (len_t i = 0; i < N_hot; i++) {
        real_t Delta;

        if (PhiRefLow[i] == 0)
            Delta = abs(PhiKinLow[i]);
        else
            Delta = abs((PhiRefLow[i]-PhiKinLow[i]) / PhiRefLow[i]);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "KineticKinetic deviates from reference flux on lower grid at index "
                LEN_T_PRINTF_FMT " with %s =/= 0. Delta = %.12e",
                i, coeffName.c_str(), Delta
            );
            success = false;
            break;
        }
    }

    // Compare fluxes (upper grid)
    for (len_t i = 0; i < N_re; i++) {
        real_t Delta;

        if (PhiRefUpp[N_hot+i] == 0)
            Delta = abs(PhiKinUpp[N_hot+i]);
        else
            Delta = abs((PhiRefUpp[N_hot+i]-PhiKinUpp[N_hot+i]) / PhiRefUpp[N_hot+i]);

        if (Delta > TOLERANCE) {
            this->PrintError(
                "KineticKinetic deviates from reference flux on upper grid at index "
                LEN_T_PRINTF_FMT " with %s =/= 0. Delta = %.12e",
                i, coeffName.c_str(), Delta
            );
            success = false;
            break;
        }
    }

    delete [] f_hot;

    delete [] PhiRefUpp;
    delete [] PhiRefLow;
    delete [] PhiKinUpp;
    delete [] PhiKinLow;

    delete refLowMat;
    delete crossLowMat;
    delete crossUpper;
    delete crossLower;
    delete refUpp;
    delete refLow;

    return success;
}*/

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
        //f_hot[i] = 1.0 + i;
		f_hot[i] = 1.0;
    for (len_t i = 0; i < N_re; i++)
        //f_re[i]  = f_hot[N_hot-1] - i*(real_t(N_hot)/real_t(N_re));
        //f_re[i] = 1.0 + i*i;
        //f_re[i] = N_hot + i;
		f_re[i] = 1.0;
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
    /*densMat->View(DREAM::FVM::Matrix::BINARY_MATLAB, "petsc_dens");*/

    // Compare lower/upper fluxes element-by-element
    //   PhiL * VpL * dxiL = sum[ PhiU * VpU * dxiBar ]
    //real_t *PhiUpper_conv = ConvertFlux(PhiUpper, runawayGrid, hottailGrid);
    real_t *PhiLower_conv = ConvertFlux(PhiLower, hottailGrid, runawayGrid);
	//if ((N_hot < N_re && N_hot%N_re == 0) || (N_hot > N_re && N_re%N_hot == 0)) {
		for (len_t i = 0; i < N_re; i++) {
			real_t Delta;

			if (PhiUpper[i] == 0)
				Delta = abs(PhiLower_conv[i]);
			else
				Delta = abs(PhiLower_conv[i] / PhiUpper[i] + 1.0);

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
	//}

    // Integrate Phi_hot and Phi_RE over momentum (p and xi)
    real_t *lowerI = hottailGrid->IntegralMomentum(PhiLower);
    real_t *upperI = runawayGrid->IntegralMomentum(PhiUpper);
	real_t *lowerI_conv = runawayGrid->IntegralMomentum(PhiLower_conv);

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
		if (upperI[ir] == 0)
			Delta = abs(lowerI_conv[ir]);
		else
			Delta = abs(lowerI_conv[ir] / upperI[ir] + 1);

		if (Delta > TOLERANCE) {
			this->PrintError(
				"hot -> re [conv]: deviation at ir = " LEN_T_PRINTF_FMT " "
				"with %s =/= 0. Delta = %e.",
				ir, coeffName.c_str(), Delta
			);
            break;
		}
    }

    delete [] lowerI;
    delete [] upperI;

    delete [] PhiLower_conv;
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
            *dxi1  = mg1->GetDp2(),
            *dxi2  = mg2->GetDp2(),
            *xi1_f = mg1->GetP2_f(),
            *xi2_f = mg2->GetP2_f();

        /*if (ir == 0) {
            SFile *sf = SFile::Create("grid.mat", SFILE_MODE_WRITE);

            sf->WriteList("Vp_re", Vp1, np1*nxi1);
            sf->WriteList("Vp_hot", Vp2, np2*nxi2);
            sf->WriteList("dp_re", dp1, np1);
            sf->WriteList("dp_hot", dp2, np2);
            sf->WriteList("dxi_re", dxi1, nxi1);
            sf->WriteList("dxi_hot", dxi2, nxi2);
            sf->WriteList("xi_re_f", xi1_f, nxi1+1);
            sf->WriteList("xi_hot_f", xi2_f, nxi2+1);

            sf->Close();
        }*/

        for (len_t j = 0; j < nxi2; j++) {
            len_t idx2   = j*np2;

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
                len_t idx1 = J*np1 + np1-1;
                real_t dxiBar = min(xi1_f[J+1], xi2_f[j+1]) - max(xi1_f[J], xi2_f[j]);

                s += Phi1[offset1+idx1] * Vp1[idx1] * dxiBar * dp1[0]
                    ;//* 2*dxi2[j]/dxiBar / (1.0+dxi2[j]/dxiBar);

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
 * Compare PXiExternalKineticKinetic to the implementations of
 * AdvectionTerm and DiffusionTerm.
 */
bool PXiExternalKineticKinetic::CompareToAdvectionDiffusionTerm() {
    return
        // Test both upwind and downwind interpolation for the
        // advection term...
        CompareToAdvectionDiffusionTerm_inner(5.4) &&
        CompareToAdvectionDiffusionTerm_inner(-5.4);
}
bool PXiExternalKineticKinetic::CompareToAdvectionDiffusionTerm_inner(const real_t coeff) {
    const len_t nr = 10, np = 4, nxi = 30;
    const len_t id_f_hot = 0, id_f_re = 1;
    bool success = true;

    DREAM::FVM::Grid *hottailGrid = this->InitializeGridRCylPXi(nr, np, nxi);
    DREAM::FVM::Grid *runawayGrid = this->InitializeGridRCylPXi(
        nr, np, nxi, hottailGrid->GetRadialGrid()->GetBmin(0),
        hottailGrid->GetMomentumGrid(0)->GetP1_f(np),  // pmin
        hottailGrid->GetMomentumGrid(0)->GetP1_f(np)*2 // pmax
    );
    DREAM::FVM::Grid *fullGrid = this->InitializeGridRCylPXi(
        nr, 2*np, nxi, hottailGrid->GetRadialGrid()->GetBmin(0),
        0, runawayGrid->GetMomentumGrid(0)->GetP1_f(np)
    );

    // Test only advection term
    DREAM::FVM::Operator *eqnBC = new DREAM::FVM::Operator(hottailGrid);
    eqnBC->AddTerm(new GeneralAdvectionTerm(hottailGrid, coeff));

    DREAM::FVM::BC::PXiExternalKineticKinetic *eqnBChot =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            hottailGrid, hottailGrid, runawayGrid,
            eqnBC, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
        );
    DREAM::FVM::BC::PXiExternalKineticKinetic *eqnBCre =
        new DREAM::FVM::BC::PXiExternalKineticKinetic(
            runawayGrid, hottailGrid, runawayGrid,
            eqnBC, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
        );

    DREAM::FVM::Operator *eqnFull = new DREAM::FVM::Operator(fullGrid);
    eqnFull->AddTerm(new GeneralAdvectionTerm(fullGrid, coeff));
    eqnFull->SetAdvectionInterpolationMethod(
        DREAM::FVM::AdvectionInterpolationCoefficient::AD_INTERP_UPWIND,
        DREAM::FVM::FLUXGRIDTYPE_P1, id_f_hot, 1.0      // 1.0 = flux limiter damping
    );

    eqnBC->RebuildTerms(101, np, nullptr);   // 101 = F1 non-zero at ip = np
    eqnFull->RebuildTerms(101, np, nullptr);

    success = success && CheckAdvectionDiffusion(
        eqnBChot, eqnBCre, eqnFull, (coeff > 0 ? "+Fp" : "-Fp"),
        hottailGrid, runawayGrid, fullGrid
    );

    delete eqnFull;
    delete eqnBCre;
    delete eqnBChot;
    delete eqnBC;

    // Test only diffusion term
    eqnBC = new DREAM::FVM::Operator(hottailGrid);
    eqnBC->AddTerm(new GeneralDiffusionTerm(hottailGrid, coeff));

    eqnBChot = new DREAM::FVM::BC::PXiExternalKineticKinetic(
            hottailGrid, hottailGrid, runawayGrid,
            eqnBC, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_LOWER
        );
    eqnBCre = new DREAM::FVM::BC::PXiExternalKineticKinetic(
            runawayGrid, hottailGrid, runawayGrid,
            eqnBC, id_f_hot, id_f_re, DREAM::FVM::BC::PXiExternalKineticKinetic::TYPE_UPPER
        );

    eqnFull = new DREAM::FVM::Operator(fullGrid);
    eqnFull->AddTerm(new GeneralDiffusionTerm(fullGrid, coeff));

    eqnBC->RebuildTerms(101, np, nullptr);   // 101 = D11 non-zero at ip = np
    eqnFull->RebuildTerms(101, np, nullptr);

    success = success && CheckAdvectionDiffusion(
        eqnBChot, eqnBCre, eqnFull, (coeff > 0 ? "+Dpp" : "-Dpp"),
        hottailGrid, runawayGrid, fullGrid
    );

    delete eqnFull;
    delete eqnBCre;
    delete eqnBChot;
    delete eqnBC;

    return success;
}

/**
 */
bool PXiExternalKineticKinetic::CheckAdvectionDiffusion(
    DREAM::FVM::BC::PXiExternalKineticKinetic *eqnBChot,
    DREAM::FVM::BC::PXiExternalKineticKinetic *eqnBCre,
    DREAM::FVM::Operator *eqnFull,
    const std::string& coeffName, DREAM::FVM::Grid *hottailGrid,
    DREAM::FVM::Grid *runawayGrid, DREAM::FVM::Grid *fullGrid
) {
    bool success = true;
    const len_t
        N_hot  = hottailGrid->GetNCells(),
        N_re   = runawayGrid->GetNCells(),
        N_full = fullGrid->GetNCells();

    DREAM::FVM::BlockMatrix *matDbl  = new DREAM::FVM::BlockMatrix();
    DREAM::FVM::BlockMatrix *matFull = new DREAM::FVM::BlockMatrix();

    matDbl->CreateSubEquation(N_hot, 5, 0);
    matDbl->CreateSubEquation(N_re,  5, 1);
    matDbl->ConstructSystem();

    matFull->CreateSubEquation(N_full, 5, 0);
    matFull->ConstructSystem();

    eqnBChot->AddToMatrixElements(matDbl, nullptr);
    matDbl->SetOffset(N_hot, N_hot);
    eqnBCre->AddToMatrixElements(matDbl, nullptr);
    matDbl->SetOffset(0, 0);

    eqnFull->SetMatrixElements(matFull, nullptr);

    matDbl->Assemble();
    matFull->Assemble();

    // Compare matrix elements
    /*PetscScalar *rowBC   = new PetscScalar[N_hot+N_re];
    PetscScalar *rowFull = new PetscScalar[N_full];
    auto getrow = [](DREAM::FVM::Matrix *mat, const PetscScalar *vals, const len_t row) {
            MatGetRow(mat->mat(), row, NULL, NULL, &vals);
            return vals;
        };*/
    auto getel = [](DREAM::FVM::Matrix *mat, PetscInt irow, PetscInt icol) {
        PetscScalar v[1];
        MatGetValues(mat->mat(), 1, &irow, 1, &icol, v);

        return v[0];
    };

    // Compare elements
    // Every row should contain exactly the same elements, just
    // in slightly different places
    #define fDbl(I,J)  getel(matDbl,  offsetDbl  + j*npDbl  + (I), offsetDbl  + j*npDbl  + (J))
    //#define fRE(I,J)   getel(matDbl,  offsetDbl  + j*npDbl  + (I), offsetDbl  + j*npDbl  + (J))
    #define fFull(I,J) getel(matFull, offsetFull + j*npFull + (I), offsetFull + j*npFull + (J))
    const real_t TOLERANCE = 100.0*std::numeric_limits<real_t>::epsilon();
    const len_t nr = fullGrid->GetNr();
    len_t offsetDbl = 0, offsetFull = 0;

    for (len_t ir = 0; ir < nr && success; ir++) {
        const len_t npDbl  = hottailGrid->GetMomentumGrid(ir)->GetNp1();
        const len_t npFull = fullGrid->GetMomentumGrid(ir)->GetNp1();
        const len_t nxi    = fullGrid->GetMomentumGrid(ir)->GetNp2();

        for (len_t j = 0; j < nxi && success; j++) {
            real_t Delta, fd, ff, fsumD=0, fsumF=0;

            //getrow(matDbl, rowBC, offsetDbl+j*npDbl+npDbl-1);

            // f_hot = c*f_hot + ...
            fd = fDbl(npDbl-1, npDbl-1);
            ff = fFull(npDbl-1, npDbl-1);
            Delta = ff==0 ? fabs(fd) : fabs(1.0 - fd/ff);
            if (Delta > TOLERANCE) {
                this->PrintError(
                    "(fHot, fHot) disagrees with full implementation (%s =/= 0) "
                    "at (ir, ixi) = (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                    "Delta = %e",
                    coeffName.c_str(), ir, j, Delta
                );
                success = false;
                break;
            }
            fsumD += fabs(fd);
            fsumF += fabs(ff);

            // f_hot = c*f_re + ...
            fd = fDbl(npDbl-1, N_hot);
            ff = fFull(npDbl-1, npDbl);
            Delta = ff==0 ? fabs(fd) : fabs(1.0 - fd/ff);
            if (Delta > TOLERANCE) {
                this->PrintError(
                    "(fHot, fRE) disagrees with full implementation (%s =/= 0) "
                    "at (ir, ixi) = (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                    "Delta = %e",
                    coeffName.c_str(), ir, j, Delta
                );
                success = false;
                break;
            }
            fsumD += fabs(fd);
            fsumF += fabs(ff);

            // f_re = c*f_hot + ...
            fd = fDbl(N_hot, npDbl-1);
            ff = fFull(npDbl, npDbl-1);
            Delta = ff==0 ? fabs(fd) : fabs(1.0 - fd/ff);
            if (Delta > TOLERANCE) {
                this->PrintError(
                    "(fRE, fHot) disagrees with full implementation (%s =/= 0) "
                    "at (ir, ixi) = (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                    "Delta = %e",
                    coeffName.c_str(), ir, j, Delta
                );
                success = false;
                break;
            }
            fsumD += fabs(fd);
            fsumF += fabs(ff);

            // f_re = c*f_re + ...
            fd = fDbl(N_hot, N_hot);
            ff = fFull(npDbl, npDbl);
            Delta = ff==0 ? fabs(fd) : fabs(1.0 - fd/ff);
            if (Delta > TOLERANCE) {
                this->PrintError(
                    "(fRE, fRE) disagrees with full implementation (%s =/= 0) "
                    "at (ir, ixi) = (" LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                    "Delta = %e",
                    coeffName.c_str(), ir, j, Delta
                );
                success = false;
                break;
            }
            fsumD += fabs(fd);
            fsumF += fabs(ff);

            // Warn in case suspiciously many elements are identically zero
            // (i.e. if all elements on a row in the matrix are zero (where
            // they should not be))
            if (fsumD == 0) {
                this->PrintWarning(
                    "All elements from PXiExternalKineticKinetic are zero in row "
                    "ir = " LEN_T_PRINTF_FMT ", ixi = " LEN_T_PRINTF_FMT ".", ir, j
                );
            }
            if (fsumF == 0) {
                this->PrintWarning(
                    "All elements from advection/diffusion are zero in row "
                    "ir = " LEN_T_PRINTF_FMT ", ixi = " LEN_T_PRINTF_FMT ".", ir, j
                );
            }
        }

        offsetDbl  += nxi*npDbl;
        offsetFull += nxi*npFull;
    }
    
    return success;
}

/**
 * Evaluate the distribution function for the advection/diffusion
 * implementation comparison.
 */
void PXiExternalKineticKinetic::CheckAdvectionDiffusion_evalF(
    DREAM::FVM::Grid *grid, std::function<real_t(const real_t, const real_t)> fEval,
    real_t *f
) {
    const len_t nr = grid->GetNr();
    len_t offs = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        DREAM::FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);

        const len_t np = mg->GetNp1(), nxi = mg->GetNp2();
        const real_t *p = mg->GetP1(), *xi = mg->GetP2();
        for (len_t j = 0; j < nxi; j++) {
            for (len_t i = 0; i < np; i++) {
                f[offs + np*j + i] = fEval(p[i], xi[j]);
            }
        }
    }
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

    // Compare to the "reference implementations"
    /*if (!CompareToReference()) {
        this->PrintError("Flux does not agree with the reference implementation.");
        success = false;
    } else
        this->PrintOK("PXiExternalKineticKinetic agrees with reference implementation.");*/

    if (!CompareToAdvectionDiffusionTerm()) {
        this->PrintError("Flux does not agree with the regular advection/diffusion implementation.");
        success = false;
    } else
        this->PrintOK("PXiExternalKineticKinetic agrees with regular advection/diffusion implementation.");

    // Check that B.C. is conservative on hot-tail, runaway and fluid grids.
    // Runaway grid should be tested with different number of xi points.
    if (!CheckConsistency())
        success = false;
    else
        this->PrintOK("Boundary condition is internally consistent.");

    return success;
}

