#include "MollerDeltaAngleKernel.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>

#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAMTESTS::_DREAM;

bool MollerDeltaAngleKernel::Run(bool) {
    bool success = true;

    if (CheckMDK_DeltaInterpolationConservation())
        this->PrintOK("Interpolated delta columns preserve normalization.");
    else {
        success = false;
        this->PrintError("Test failed: interpolated delta columns do not normalize to 1.");
    }
    if (CheckMDK_KinematicEdgeCases())
        this->PrintOK("Angle kernel handles low energies, empty support and singleton ranges.");
    else {
        success = false;
        this->PrintError("Test failed: angle kernel kinematic edge cases.");
    }
    if (CheckMDK_IsotropicAndDipoleMoments())
        this->PrintOK("Isotropic input has zero drift; dipole transfer follows xi_star.");
    else {
        success = false;
        this->PrintError("Test failed: angle kernel isotropic/dipole moments.");
    }
    return success;
}

bool MollerDeltaAngleKernel::CheckMDK_IsotropicAndDipoleMoments() {
    // A cylinder isolates angular redistribution from trapped-orbit conventions.
    // Include an odd pitch grid (xi=0 cell) and unequal primary/secondary grids.
    const len_t sizes[][2] = {{30, 40}, {45, 45}, {60, 60}};
    bool success = true;
    for (const auto &size : sizes) {
        const len_t nP = size[0], nK = size[1];
        std::unique_ptr<DREAM::FVM::Grid> target(
            InitializeGridRCylPXi(1, 19, nK, 1.0, 0.0, 1.9));
        std::unique_ptr<DREAM::FVM::Grid> primary(
            InitializeGridRCylPXi(1, 19, nP, 1.0, 0.0, 1.9));
        auto *mgK = target->GetMomentumGrid(0), *mgP = primary->GetMomentumGrid(0);
        DREAM::MollerEnergyKernel energy(target.get(), primary.get(), 0.1);
        DREAM::MollerDeltaAngleKernel angle(target.get(), primary.get(), 0.1,
            80, 8);
        std::vector<real_t> isotropic(nP), dipole(nP), output(nK), outputDipole(nK);
        real_t incidentDipoleMoment = 0;
        for (len_t l=0; l<nP; ++l) {
            const real_t xi = mgP->GetXi0(0, l), dx = mgP->GetDp2(l);
            isotropic[l] = dx/2;
            dipole[l] = dx*xi;
            incidentDipoleMoment += dx*xi*xi;
        }
        // The secondary cell-centre moment is a quadrature approximation.
        // Keep its tolerance separate from exact reflection/zero-current checks.
        const real_t dxK = 2.0/nK;
        const real_t momentTolerance = 2*dxK*dxK;
        len_t active = 0;
        for (len_t i=0; i<19; ++i) for (len_t k=0; k<19; ++k) {
            if (energy.DifferentialCS(i,k) == 0) continue;
            ++active;
            std::fill(output.begin(), output.end(), 0);
            std::fill(outputDipole.begin(), outputDipole.end(), 0);
            angle.AccumulatePitch(0, i, k, isotropic.data(), 1, output.data());
            angle.AccumulatePitch(0, i, k, dipole.data(), 1, outputDipole.data());
            real_t norm = 0, drift = 0, dipoleNumber = 0, dipoleMoment = 0;
            real_t evenError = 0, oddError = 0;
            for (len_t j=0; j<nK; ++j) {
                const real_t dx = mgK->GetDp2(j), xi = mgK->GetXi0(0,j);
                norm += dx*output[j];
                drift += dx*xi*output[j];
                dipoleNumber += dx*outputDipole[j];
                dipoleMoment += dx*xi*outputDipole[j];
                evenError += dx*fabs(output[j]-output[nK-1-j]);
                oddError += dx*fabs(outputDipole[j]+outputDipole[nK-1-j]);
            }
            const real_t p1 = mgP->GetP1(k);
            const real_t lower = std::max(0.1, mgK->GetP1_f(i));
            const real_t upper = std::min(mgK->GetP1_f(i+1),
                DREAM::KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1));
            const real_t p = std::max(lower, std::min(mgK->GetP1(i), upper));
            const real_t xs = DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(p,p1);
            const real_t error = fabs(dipoleMoment/incidentDipoleMoment-xs);
            if (!std::isfinite(norm) || !std::isfinite(error)
                || fabs(norm-1)>5e-13 || fabs(drift)>5e-13
                || fabs(dipoleNumber)>5e-13 || evenError>5e-12 || oddError>5e-12
                || error>momentTolerance) {
                this->PrintError("Angular moments failed: nP=%ld nK=%ld i=%ld k=%ld "
                    "norm=%.16g drift=%.16g dipole_error=%.16g",
                    (long)nP, (long)nK, (long)i, (long)k, norm, drift, error);
                success = false;
            }
        }
        success &= active>0;
    }
    return success;
}

bool MollerDeltaAngleKernel::CheckMDK_KinematicEdgeCases() {
    struct Case { len_t npK; real_t maxK; len_t npP; real_t maxP; real_t cutoff; bool empty; };
    const Case cases[] = {
        {10, 1.0, 10, 1.0, 0.1, false},  // all primaries below sqrt(3)
        {10, 1.0, 10, 1.0, 2.0, true},   // no kinematically accessible cells
        {1, 4.0, 1, 4.0, 0.2, false},    // single scattering angle
        {1, 4.0, 10, 1.0, 0.1, false},   // target centre above primary momentum
        {19, 1.9, 19, 1.9, 0.1, false}   // old range silently clamps valid angles
    };
    bool success = true;
    for (const auto &c : cases) {
        std::unique_ptr<DREAM::FVM::Grid> target(
            InitializeGridRCylPXi(1, c.npK, 40, 1.0, 0.0, c.maxK));
        std::unique_ptr<DREAM::FVM::Grid> primary(
            InitializeGridRCylPXi(1, c.npP, 20, 1.0, 0.0, c.maxP));
        DREAM::MollerEnergyKernel energy(target.get(), primary.get(), c.cutoff);
        DREAM::MollerDeltaAngleKernel angle(target.get(), primary.get(), c.cutoff,
            80, 8);
        auto *mgK = target->GetMomentumGrid(0);
        auto *mgP = primary->GetMomentumGrid(0);
        std::vector<real_t> weights(20, 0), pitch(40, 0), direct(40, 0);
        weights[19] = 1;
        len_t active = 0;
        for (len_t i=0; i<c.npK; ++i) for (len_t k=0; k<c.npP; ++k) {
            if (energy.DifferentialCS(i,k) == 0) continue;
            ++active;
            std::fill(pitch.begin(), pitch.end(), 0);
            angle.AccumulatePitch(0, i, k, weights.data(), 1, pitch.data());
            const real_t p1 = mgP->GetP1(k);
            const real_t plo = std::max(c.cutoff, mgK->GetP1_f(i));
            const real_t phi = std::min(mgK->GetP1_f(i+1),
                DREAM::KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1));
            const real_t p = std::max(plo, std::min(mgK->GetP1(i), phi));
            const real_t xs = DREAM::KnockOnUtilities::Kinematics::EvaluateXiStar(p, p1);
            DREAM::KnockOnUtilities::SetDeltaMatrixColumnOnGrid(0, xs, 19,
                target.get(), primary.get(), direct.data(), 8);
            real_t norm = 0, error = 0;
            for (len_t j=0; j<40; ++j) {
                success &= std::isfinite(pitch[j]);
                norm += mgK->GetDp2(j)*pitch[j];
                error += mgK->GetDp2(j)*fabs(pitch[j]-direct[j]);
            }
            success &= std::isfinite(norm) && fabs(norm-1) < 1e-12;
            // Allow interpolation error at the sharp edges of angular support;
            // the old range gives order-unity errors in this comparison.
            success &= std::isfinite(error) && error < 0.13;
        }
        success &= c.empty ? (active == 0) : (active > 0);
        if (c.empty) {
            angle.AccumulatePitch(0, 0, 0, weights.data(), 1, pitch.data());
            for (auto x : pitch) success &= std::isfinite(x) && x == 0;
        }
    }
    return success;
}

bool MollerDeltaAngleKernel::CheckMDK_DeltaInterpolationConservation() {
    const real_t tol = 1e-14;

    len_t nr = 2;
    len_t npK = 5;
    len_t nxiK = 40;
    len_t npP = 7;
    len_t nxiP = 30;

    len_t ntheta_interp = 50;
    len_t nrProfiles = 8;
    real_t pMin = 0;
    real_t pMax = 3;

    real_t B0 = 1.0;
    auto *gridF = InitializeFluidGrid(nr, B0);
    auto *gridK = InitializeGridGeneralRPXi(nr, npK, nxiK, ntheta_interp, nrProfiles, pMin, pMax);
    auto *gridP = InitializeGridGeneralRPXi(nr, npP, nxiP, ntheta_interp, nrProfiles, pMin, pMax);

    const real_t pCutoff = gridK->GetMomentumGrid(0)->GetP1(1);
    constexpr len_t n_xi_stars_tabulate = 80;
    constexpr len_t poloidal_order = 8;

    DREAM::MollerDeltaAngleKernel K(
        gridK, gridP, pCutoff, n_xi_stars_tabulate, poloidal_order
    );
    K.GridRebuilt();

    bool success = true;

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mgK = gridK->GetMomentumGrid(ir);
        auto *mgP = gridP->GetMomentumGrid(ir);

        const len_t NpK = mgK->GetNp1();
        const len_t NxiK = mgK->GetNp2();
        const len_t NpP = mgP->GetNp1();
        const len_t NxiP = mgP->GetNp2();

        std::vector<real_t> W_l(NxiP, 0.0);
        std::vector<real_t> outPitch(NxiK, 0.0);

        for (len_t i = 0; i < NpK; i++) {
            for (len_t k = 0; k < NpP; k++) {
                for (len_t l = 0; l < NxiP; l++) {
                    std::fill(W_l.begin(), W_l.end(), 0.0);
                    W_l[l] = 1.0;

                    std::fill(outPitch.begin(), outPitch.end(), 0.0);
                    K.AccumulatePitch(ir, i, k, W_l.data(), /*Sik=*/1.0, outPitch.data());

                    real_t sum = 0.0;
                    for (len_t j = 0; j < NxiK; j++) sum += mgK->GetDp2(j) * outPitch[j];

                    const bool primaryVoid = (gridP->GetVpOverP2AtZero(ir)[l] == 0);
                    if (primaryVoid) {
                        if (fabs(sum) != 0) {
                            this->PrintError(
                                "Delta normalization non-zero in void primary cell at ir=%ld i=%ld "
                                "k=%ld l=%ld: sum=%.16g",
                                ir, i, k, l, sum
                            );
                            success = false;
                        }
                    } else {
                        if (!std::isfinite(sum) || fabs(sum - 1.0) > tol) {
                            this->PrintError(
                                "Delta normalization failed at ir=%ld i=%ld k=%ld l=%ld: sum=%.16g",
                                ir, i, k, l, sum
                            );
                            success = false;
                        }
                    }
                }
            }
        }
    }

    delete gridF;
    delete gridK;
    delete gridP;

    return success;
}
