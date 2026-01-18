/**
 * Tests geometry calculations in DREAM::FVM::RadialGridGenerator using
 * different radial grid generators (e.g. AnalyticB and NumericB).
 */
#include "RadialGridGenerator.hpp"

#include <cmath>
#include <string>

#include "FVM/Grid/NumericBRadialGridGenerator.hpp"
#include "UnitTest.hpp"

using namespace DREAMTESTS::FVM;

/**
 * Compare two floating-point values with an atol and rtol using the same
 * conventions as, for example, in numpy and scipy.
 */
static bool nearlyEqual(real_t a, real_t b, real_t atol, real_t rtol) {
    const real_t diff = std::fabs(a - b);
    const real_t scale = std::max(std::fabs(a), std::fabs(b));
    return diff <= atol + rtol * scale;
}

/**
 * Print a detailed comparison of two floating-point values, including
 * absolute/relative differences and the effective tolerance.
 *
 * Used to provide diagnostics when a test comparison fails.
 */
void RadialGridGenerator::printComparison(
    real_t val1, real_t val2, const std::string &name1, const std::string &name2, real_t atol,
    real_t rtol, len_t n_indent
) {
    const real_t adiff = std::fabs(val1 - val2);
    const real_t scale = std::max(std::fabs(val1), std::fabs(val2));

    std::string indent = std::string(n_indent, ' ');
    this->PrintError(indent + name1 + " = %.8g\n", val1);
    this->PrintError(indent + name2 + " = %.8g\n", val2);
    this->PrintError(indent + "|diff|    = %.8g\n", adiff);

    if (scale > 0) {
        this->PrintError(indent + "rel diff  = %.8g\n", adiff / scale);
    }
    this->PrintError(
        indent + "tol      = atol + rtol*scale = %.8g + %.8g*%.8g = %.8g\n", atol, rtol, scale,
        atol + rtol * scale
    );
}

/**
 * At a single (r, theta) location, verify that individual geometry calculations
 * (BAtTheta, JacobianAtTheta, etc.) agree with the combined
 * GeometricQuantitiesAtTheta() evaluation.
 */
bool RadialGridGenerator::TestGeometricQuantitiesAtTheta(
    len_t ir, real_t theta, DREAM::FVM::RadialGrid *rg, real_t atol, real_t rtol
) {
    auto *fsa = rg->GetFluxSurfaceAverager();

    bool success = true;
    real_t B1, Jacobian1, ROverR01, NablaR21;
    fsa->GeometricQuantitiesAtTheta(
        ir, theta, B1, Jacobian1, ROverR01, NablaR21, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
    );
    real_t B2 = fsa->BAtTheta(ir, theta, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t Jacobian2 = fsa->JacobianAtTheta(ir, theta, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t ROverR02 = fsa->ROverR0AtTheta(ir, theta, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t NablaR22 = fsa->NablaR2AtTheta(ir, theta, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION);
    if (!nearlyEqual(B1, B2, atol, rtol)) {
        this->PrintError("BAtTheta does not agree with GeometricQuantitiesAtTheta.\n");
        printComparison(B1, B2, "GeometricQuantitiesAtTheta", "BAtTheta", atol, rtol);
        success = false;
    }
    if (!nearlyEqual(Jacobian1, Jacobian2, atol, rtol)) {
        this->PrintError("JacobianAtTheta does not agree with GeometricQuantitiesAtTheta.\n");
        printComparison(
            Jacobian1, Jacobian2, "GeometricQuantitiesAtTheta", "JacobianAtTheta", atol, rtol
        );
        success = false;
    }
    if (!nearlyEqual(ROverR01, ROverR02, atol, rtol)) {
        this->PrintError("ROverR0AtTheta does not agree with GeometricQuantitiesAtTheta.\n");
        printComparison(
            ROverR01, ROverR02, "GeometricQuantitiesAtTheta", "ROverR0AtTheta", atol, rtol
        );
        success = false;
    }
    if (!nearlyEqual(NablaR21, NablaR22, atol, rtol)) {
        this->PrintError("NablaR2AtTheta does not agree with GeometricQuantitiesAtTheta.\n");
        printComparison(
            NablaR21, NablaR22, "GeometricQuantitiesAtTheta", "NablaR2AtTheta", atol, rtol
        );
        success = false;
    }
    return success;
}

/**
 * Verify internal consistency of geometry calculations in the NumericB grid:
 * individual methods (BAtTheta, JacobianAtTheta, etc.) must agree with
 * GeometricQuantitiesAtTheta().
 */
bool RadialGridGenerator::CheckNumericBGeometryCalculations() {
    real_t atol = 5 * std::numeric_limits<real_t>::epsilon();
    real_t rtol = 5 * std::numeric_limits<real_t>::epsilon();
    len_t nr = 4;
    len_t ntheta_interp = 50;
    auto *nbrgg = this->InitializeNumericBRadialGridGenerator(nr, ntheta_interp);
    auto *rg = new DREAM::FVM::RadialGrid(nbrgg, 0);
    rg->Rebuild(0);
    bool success = true;

    len_t ntheta = 10;
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t n = 0; n < ntheta; n++) {
            real_t theta = 2 * M_PI * n / (ntheta - 1);
            success &= TestGeometricQuantitiesAtTheta(ir, theta, rg, atol, rtol);
        }
    }
    delete rg;
    return success;
}

/**
 * Verify internal consistency of geometry calculations in the AnalyticB grid:
 * individual methods (BAtTheta, JacobianAtTheta, etc.) must agree with
 * GeometricQuantitiesAtTheta().
 */
bool RadialGridGenerator::CheckAnalyticBGeometryCalculations() {
    // Tolerance slightly above machine epsilon to account for floating-point roundoff
    real_t atol = 5 * std::numeric_limits<real_t>::epsilon();
    real_t rtol = 5 * std::numeric_limits<real_t>::epsilon();
    len_t nr = 4;
    len_t ntheta_interp = 50;
    len_t nrProfiles = 20;
    auto *abrgg = this->InitializeAnalyticBRadialGridGenerator(nr, nrProfiles, ntheta_interp);
    auto *rg = new DREAM::FVM::RadialGrid(abrgg, 0);
    rg->Rebuild(0);
    bool success = true;

    len_t ntheta = 10;
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t n = 0; n < ntheta; n++) {
            real_t theta = 2 * M_PI * n / (ntheta - 1);
            success &= TestGeometricQuantitiesAtTheta(ir, theta, rg, atol, rtol);
        }
    }
    delete rg;
    return success;
}

/**
 * Creates NumericB and AnalyticB radial grids that describe the same
 * underlying magnetic field, and verifies that geometric calculations
 * (Jacobian, B, ...) agree within small error at all radii and poloidal locations.
 */
bool RadialGridGenerator::CheckNumericBAgreesWithAnalyticB() {
    // A relatively loose tolerance is used, so that we can
    // get away with using a smaller dataset. Higher-resolution
    // magnetic field data yields smaller errors.
    real_t atol = 1e-6;
    real_t rtol = 5e-4;

    len_t nr = 6;
    len_t ntheta_interp = 200;

    // the DeltaPrime value of the example grid
    real_t DeltaPrime = 0.3;

    auto *nbrgg = this->InitializeNumericBRadialGridGenerator(nr, ntheta_interp);
    auto *rg_n = new DREAM::FVM::RadialGrid(nbrgg, 0);
    rg_n->Rebuild(0);
    auto *fsa_n = rg_n->GetFluxSurfaceAverager();

    len_t nrProfiles = 20;
    auto *abrgg = this->InitializeAnalyticBRadialGridGenerator(nr, nrProfiles, ntheta_interp);
    auto *rg_a = new DREAM::FVM::RadialGrid(abrgg, 0);
    rg_a->Rebuild(0);
    auto *fsa_a = rg_a->GetFluxSurfaceAverager();

    bool success = true;
    len_t ntheta = 10;
    for (len_t ir = 0; ir < nr; ir++) {
        for (len_t n = 0; n < ntheta; n++) {
            real_t theta = 2 * M_PI * n / (ntheta - 1);
            real_t B_a, Jacobian_a, ROverR0_a, NablaR2_a;
            fsa_a->GeometricQuantitiesAtTheta(
                ir, theta, B_a, Jacobian_a, ROverR0_a, NablaR2_a,
                DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
            );
            real_t B_n, Jacobian_n, ROverR0_n, NablaR2_n;
            fsa_n->GeometricQuantitiesAtTheta(
                ir, theta, B_n, Jacobian_n, ROverR0_n, NablaR2_n,
                DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
            );

            // Compensate for differences in radial coordinate definitions:
            //   - Numeric LUKE equilibrium uses r = R - Rp on the outer midplane
            //   - Analytic geometry uses R = Rp + Delta(r) + r
            //
            // This introduces a radially varying Shafranov shift, meaning the same (R, Z)
            // location corresponds to different minor radii r in the two coordinate systems.
            real_t B_pred = B_n;
            real_t ROverR0_pred = ROverR0_n;
            real_t Jacobian_pred = Jacobian_n * (1 + DeltaPrime);
            real_t Nabla2_pred = NablaR2_n / ((1 + DeltaPrime) * (1 + DeltaPrime));
            if (!nearlyEqual(B_pred, B_a, atol, rtol)) {
                this->PrintError(
                    "Evaluation of B does not agree between NumericB and AnalyticB at ir=%ld "
                    "(r=%.4g), theta=%.4g.\n",
                    ir, rg_n->GetR(ir), theta
                );
                printComparison(B_pred, B_a, "NumericB", "AnalyticB", atol, rtol);
                success = false;
            }
            if (!nearlyEqual(Jacobian_pred, Jacobian_a, atol, rtol)) {
                this->PrintError(
                    "Evaluation of Jacobian does not agree between NumericB and AnalyticB at "
                    "ir=%ld (r=%.4g), theta=%.4g.\n",
                    ir, rg_n->GetR(ir), theta
                );
                printComparison(Jacobian_pred, Jacobian_a, "NumericB", "AnalyticB", atol, rtol);
                success = false;
            }
            if (!nearlyEqual(ROverR0_pred, ROverR0_a, atol, rtol)) {
                this->PrintError(
                    "Evaluation of ROverR0 does not agree between NumericB and AnalyticB at ir=%ld "
                    "(r=%.4g), theta=%.4g.\n",
                    ir, rg_n->GetR(ir), theta
                );
                printComparison(ROverR0_pred, ROverR0_a, "NumericB", "AnalyticB", atol, rtol);
                success = false;
            }
            if (!nearlyEqual(Nabla2_pred, NablaR2_a, atol, rtol)) {
                this->PrintError(
                    "Evaluation of NablaR2 does not agree between NumericB and AnalyticB at ir=%ld "
                    "(r=%.4g), theta=%.4g.\n",
                    ir, rg_n->GetR(ir), theta
                );
                printComparison(Nabla2_pred, NablaR2_a, "NumericB", "AnalyticB", atol, rtol);
                success = false;
            }
        }
    }
    delete rg_a;
    delete rg_n;
    return success;
}

/**
 * Run all RadialGridGenerator tests.
 *
 * Returns true if all tests pass, false otherwise.
 */
bool RadialGridGenerator::Run(bool) {
    bool success = true;

    if (CheckAnalyticBGeometryCalculations())
        this->PrintOK("Successfully checked geometry calculations in the AnalyticB radial grid.");
    else {
        success = false;
        this->PrintError("Invalid geometry calculations in the AnalyticB radial grid.");
    }
    if (CheckNumericBGeometryCalculations())
        this->PrintOK("Successfully checked geometry calculations in the NumericB radial grid.");
    else {
        success = false;
        this->PrintError("Invalid geometry calculations in the NumericB radial grid.");
    }
    if (CheckNumericBAgreesWithAnalyticB())
        this->PrintOK(
            "NumericB and AnalyticB give the same geometric quantities in the same magnetic field."
        );
    else {
        success = false;
        this->PrintError(
            "NumericB and AnalyticB give different geometric quantities in the same magnetic field."
        );
    }

    return success;
}
