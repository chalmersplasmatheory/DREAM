/**
 * Implements tests of geometry calculations in DREAM::FVM::RadialGridGenerator
 * with different radial grid generators, such as Cylindrical, AnalyticB and Numerical.
 */
#include "RadialGridGenerator.hpp"
#include "FVM/Grid/NumericBRadialGridGenerator.hpp"
#include "UnitTest.hpp"

using namespace DREAMTESTS::FVM;

namespace {
bool TestGeometricQuantitiesAtTheta(
    len_t ir, real_t theta, DREAM::FVM::RadialGrid *rg, real_t test_tol
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
    if (fabs(B1 - B2) > test_tol) {
        printf("Evaluation of B does not agree between BAtTheta and GeometricQuantitiesAtTheta.\n");
        printf("B (GeometricQuantitiesAtTheta): %.4g\n", B1);
        printf("B (BAtTheta): %.4g\n", B2);
        success = false;
    }
    if (fabs(Jacobian1 - Jacobian2) > test_tol) {
        printf("Evaluation of Jacobian does not agree between JacobianAtTheta and "
               "GeometricQuantitiesAtTheta.\n");
        printf("Jacobian (GeometricQuantitiesAtTheta): %.4g\n", Jacobian1);
        printf("Jacobian (JacobianAtTheta): %.4g\n", Jacobian2);
        success = false;
    }
    if (fabs(ROverR01 - ROverR02) > test_tol) {
        printf("Evaluation of ROverR0 does not agree between ROverR0AtTheta and "
               "GeometricQuantitiesAtTheta.\n");
        printf("ROverR0 (GeometricQuantitiesAtTheta): %.4g\n", ROverR01);
        printf("ROverR0 (ROverR0AtTheta): %.4g\n", ROverR02);
        success = false;
    }
    if (fabs(NablaR21 - NablaR22) > test_tol) {
        printf("Evaluation of NablaR2 does not agree between NablaR2AtTheta and "
               "GeometricQuantitiesAtTheta.\n");
        printf("NablaR2 (GeometricQuantitiesAtTheta): %.4g\n", NablaR21);
        printf("NablaR2 (NablaR2AtTheta): %.4g\n", NablaR22);
        success = false;
    }
    return success;
}
} // namespace

bool RadialGridGenerator::CheckNumericBGeometryCalculations() {
    real_t test_tol = 1e-10;
    len_t nr = 4;
    real_t ntheta_interp = 50;
    auto *nbrgg = this->InitializeNumericBRadialGridGenerator(nr, ntheta_interp);
    auto *rg = new DREAM::FVM::RadialGrid(nbrgg, 0);
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
    len_t ntheta_interp = 50;
    len_t nrProfiles = 20;
    auto *abrgg = this->InitializeAnalyticBRadialGridGenerator(nr, nrProfiles, ntheta_interp);
    auto *rg = new DREAM::FVM::RadialGrid(abrgg, 0);
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

bool RadialGridGenerator::CheckNumericBAgreesWithAnalyticB(){
    // A relatively loose tolerance is used, so that we can
    // get away with using a smaller dataset. Higher-resolution
    // data yields smaller errors.
    real_t test_tol = 1e-4;

    len_t nr = 4;
    real_t ntheta_interp = 200;

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
        real_t theta = 0;
        for (len_t n = 0; n < ntheta; n++) {
            real_t B_a, Jacobian_a, ROverR0_a, NablaR2_a;
            fsa_a->GeometricQuantitiesAtTheta(
                ir, theta, B_a, Jacobian_a, ROverR0_a, NablaR2_a, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
            );
            real_t B_n, Jacobian_n, ROverR0_n, NablaR2_n;
            fsa_n->GeometricQuantitiesAtTheta(
                ir, theta, B_n, Jacobian_n, ROverR0_n, NablaR2_n, DREAM::FVM::FLUXGRIDTYPE_DISTRIBUTION
            );

            // Now compensate for the fact that the numeric LUKE equilibrium uses: 
            //    r = R - Rp
            // on the outer midplane as radial parameter, whereas in the analytic geometry:
            //    R = Rp + Delta(r) + r
            // on the outer midplane, e.g. there is a radially varying shafranov shift difference.
            real_t B_pred = B_n;
            real_t ROverR0_pred = ROverR0_n;
            real_t Jacobian_pred = Jacobian_n * (1 + DeltaPrime);
            real_t Nabla2_pred = NablaR2_n / ((1 + DeltaPrime) * (1 + DeltaPrime));
            if( fabs(B_pred - B_a) > test_tol ){
                printf("Evaluation of B does not agree between NumericB and AnalyticB at ir=%ld (r=%.4g), theta=%.4g.\n", ir, rg_n->GetR(ir), theta);
                printf("B (Numeric): %.4g\n", B_pred);
                printf("B (Analytic): %.4g\n", B_a);
                success = false;
            }
            if( fabs(Jacobian_pred - Jacobian_a) > test_tol ){
                printf("Evaluation of Jacobian does not agree between NumericB and AnalyticB at ir=%ld (r=%.4g), theta=%.4g.\n", ir, rg_n->GetR(ir), theta);
                printf("Jacobian (Numeric): %.4g\n", Jacobian_pred);
                printf("Jacobian (Analytic): %.4g\n", Jacobian_a);
                success = false;
            }
            if( fabs(ROverR0_pred - ROverR0_a) > test_tol ){
                printf("Evaluation of ROverR0 does not agree between NumericB and AnalyticB at ir=%ld (r=%.4g), theta=%.4g.\n", ir, rg_n->GetR(ir), theta);
                printf("ROverR0 (Numeric): %.4g\n", ROverR0_pred);
                printf("ROverR0 (Analytic): %.4g\n", ROverR0_a);
                success = false;
            }
            if( fabs(Nabla2_pred - NablaR2_a) > test_tol ){
                printf("Evaluation of NablaR2 does not agree between NumericB and AnalyticB at ir=%ld (r=%.4g), theta=%.4g.\n", ir, rg_n->GetR(ir), theta);
                printf("NablaR2 (Numeric): %.4g\n", Nabla2_pred);
                printf("NablaR2 (Analytic): %.4g\n", NablaR2_a);
                success = false;
            }
            theta += 2 * M_PI / (ntheta - 1);
        }
    }
    delete rg_a;
    delete rg_n;
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
    if (CheckNumericBAgreesWithAnalyticB())
        this->PrintOK("NumericB and AnalyticB give the same geometric quantities in the same magnetic field.");
    else {
        success = false;
        this->PrintError("NumericB and AnalyticB give different geometric quantities in the same magnetic field.");
    }

    
    return success;
}
