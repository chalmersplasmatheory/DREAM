/**
 * Implements tests of geometry calculations in DREAM::FVM::RadialGridGenerator
 * with different radial grid generators, such as Cylindrical, AnalyticB and Numerical.
 */
#include "RadialGridGenerator.hpp"
#include "UnitTest.hpp"
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"

using namespace DREAMTESTS::FVM;

namespace {
    bool TestGeometricQuantitiesAtTheta(len_t ir, real_t theta, DREAM::FVM::RadialGrid *rg, real_t test_tol){
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
        if (fabs(B1 - B2) > test_tol){
            printf("Evaluation of B does not agree between BAtTheta and GeometricQuantitiesAtTheta.\n");
            printf("B (GeometricQuantitiesAtTheta): %.4g\n", B1);
            printf("B (BAtTheta): %.4g\n", B2);
            success = false;
        }
        if (fabs(Jacobian1 - Jacobian2) > test_tol){
            printf("Evaluation of Jacobian does not agree between JacobianAtTheta and GeometricQuantitiesAtTheta.\n");
            printf("Jacobian (GeometricQuantitiesAtTheta): %.4g\n", Jacobian1);
            printf("Jacobian (JacobianAtTheta): %.4g\n", Jacobian2);
            success = false;
        }
        if (fabs(ROverR01 - ROverR02) > test_tol){
            printf("Evaluation of ROverR0 does not agree between ROverR0AtTheta and GeometricQuantitiesAtTheta.\n");
            printf("ROverR0 (GeometricQuantitiesAtTheta): %.4g\n", ROverR01);
            printf("ROverR0 (ROverR0AtTheta): %.4g\n", ROverR02);
            success = false;
        }
        if (fabs(NablaR21 - NablaR22) > test_tol){
            printf("Evaluation of NablaR2 does not agree between NablaR2AtTheta and GeometricQuantitiesAtTheta.\n");
            printf("NablaR2 (GeometricQuantitiesAtTheta): %.4g\n", NablaR21);
            printf("NablaR2 (NablaR2AtTheta): %.4g\n", NablaR22);
            success = false;
        }
        return success;
    }
}

bool RadialGridGenerator::CheckNumericBGeometryCalculations() {
    real_t test_tol = 1e-10;
    len_t nr = 4;
    real_t a = 0.577;
    len_t ntheta_interp = 100;
    auto *NBrgg = new DREAM::FVM::NumericBRadialGridGenerator(
        nr, 0, a, "numericmag_physics_test_field.h5", DREAM::FVM::NumericBRadialGridGenerator::FILE_FORMAT_LUKE, ntheta_interp 
    );
    auto *rg = new DREAM::FVM::RadialGrid(NBrgg, 0);
    rg->Rebuild(0);
    bool success = true;
    len_t ntheta = 10;
    for (len_t ir = 0; ir < nr; ir++) {
        real_t theta = 0;
        for (len_t n = 0; n < ntheta; n++) {
            success = TestGeometricQuantitiesAtTheta(ir, theta, rg, test_tol);
            theta += 2 * M_PI / (ntheta - 1);
        }
    }
    delete rg;
    return success;
}


bool RadialGridGenerator::CheckAnalyticBGeometryCalculations() {
    real_t test_tol = 1e-10;
    len_t nr = 4;
    len_t ntheta_interp = 20;
    len_t nrProfiles = 20;
    auto *ABrgg = InitializeAnalyticBRadialGridGenerator(nr, nrProfiles, ntheta_interp);
    auto *rg = new DREAM::FVM::RadialGrid(ABrgg, 0);
    rg->Rebuild(0);
    bool success = true;
    len_t ntheta = 10;
    for (len_t ir = 0; ir < nr; ir++) {
        real_t theta = 0;
        for (len_t n = 0; n < ntheta; n++) {
            success = TestGeometricQuantitiesAtTheta(ir, theta, rg, test_tol);
            theta += 2 * M_PI / (ntheta - 1);
        }
    }
    delete rg;
    return success;
}

/**
 * Run all RadialGridGenerator tests.
 * Returns 'true' if all tests passed. 'false' otherwise.
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

    return success;
}
