/**
 * Tests for the 'Interpolator3D' class.
 */

#include <cmath>
#include "FVM/Interpolator3D.hpp"
#include "Interpolator3D.hpp"

using namespace DREAMTESTS::FVM;
using namespace std;


struct Interpolator3D::griddata *Interpolator3D::GenerateData(
    struct gridlimits& limits,
    const len_t nx1, const len_t nx2, const len_t nx3,
    function<real_t(real_t, real_t, real_t)>& func
) {
    struct griddata *gd = new struct griddata;

    gd->nx1 = nx1;
    gd->nx2 = nx2;
    gd->nx3 = nx3;

    real_t 
        *x1 = new real_t[nx1],
        *x2 = new real_t[nx2],
        *x3 = new real_t[nx3],
        *data = new real_t[nx1*nx2*nx3];

    for (len_t i = 0; i < nx1; i++)
        x1[i] = limits.x1min + (limits.x1max-limits.x1min) * i / (nx1-1);
    for (len_t i = 0; i < nx2; i++)
        x2[i] = limits.x2min + (limits.x2max-limits.x2min) * i / (nx2-1);
    for (len_t i = 0; i < nx3; i++)
        x3[i] = limits.x3min + (limits.x3max-limits.x3min) * i / (nx3-1);

    #define IDX(X1,X2,X3) (((X1)*nx2 + (X2))*nx3 + (X3))
    for (len_t i = 0; i < nx1; i++)
        for (len_t j = 0; j < nx2; j++)
            for (len_t k = 0; k < nx3; k++)
                data[IDX(i,j,k)] = func(x1[i], x2[j], x3[k]);
    #undef IDX

    gd->x1 = x1;
    gd->x2 = x2;
    gd->x3 = x3;
    gd->y  = data;

    return gd;
}


/**
 * Evaluate the given Interpolator3D on a test grid,
 * pretending to use both P/XI and PPAR/PPERP coordinates.
 */
bool Interpolator3D::EvalInterpolator3D(
    DREAM::FVM::Interpolator3D *interp,
    function<real_t(real_t, real_t, real_t)>& func
) {
    const len_t ntests = 8;
    const len_t nx1 = ntests-2, nx2 = ntests-1, nx3 = ntests;
    const real_t xtest[3][ntests] = {
        {0.09196,0.03264,0.79261,0.34556,0.54265,0.18482,0.24265,0.80674},
        {0.86540,0.14284,0.43475,0.81089,0.47123,0.08429,0.76421,0.90060},
        {0.24540,0.58885,0.71874,0.13600,0.42034,0.72153,0.62896,0.51662}
    };

    // Evaluate interpolator object
    const real_t *dataPXI = interp->Eval(
        nx1, nx2, nx3, xtest[0], xtest[1], xtest[2],
        DREAM::FVM::Interpolator3D::GRID_PXI
    );

    const real_t *dataPP  = interp->Eval(
        nx1, nx2, nx3, xtest[0], xtest[1], xtest[2],
        DREAM::FVM::Interpolator3D::GRID_PPARPPERP
    );

    // Check the interpolated values
    const real_t TOL = 1e3*std::numeric_limits<real_t>::epsilon();
    if (interp->GetGridType() == DREAM::FVM::Interpolator3D::GRID_PXI) {
        const real_t *r  = xtest[0];
        const real_t *xi = xtest[1];
        const real_t *p  = xtest[2];

        for (len_t i = 0; i < nx1; i++) {
            for (len_t j = 0; j < nx2; j++) {
                for (len_t k = 0; k < nx3; k++) {
                    len_t idx = ((i*nx2 + j)*nx3 + k);
                    real_t Delta = abs(func(r[i], xi[j], p[k]) / dataPXI[idx] - 1);

                    if (Delta > TOL) {
                        this->PrintError(
                            "3D interpolation failed for p/xi to p/xi at element ("
                            LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                            "Delta = %e",
                            i, j, k, Delta
                        );
                        return false;
                    }
                }
            }
        }

        const real_t *pperp = xtest[1];
        const real_t *ppar  = xtest[2];

        for (len_t i = 0; i < nx1; i++) {
            for (len_t j = 0; j < nx2; j++) {
                for (len_t k = 0; k < nx3; k++) {
                    len_t idx = ((i*nx2 + j)*nx3 + k);
                    real_t p  = sqrt(ppar[k]*ppar[k] + pperp[j]*pperp[j]);
                    real_t xi = ppar[k] / p;

                    real_t Delta = abs(func(r[i], xi, p) / dataPP[idx] - 1);

                    if (Delta > TOL) {
                        this->PrintError(
                            "3D interpolation failed for p/xi to ppar/pperp at element ("
                            LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                            "Delta = %e",
                            i, j, k, Delta
                        );
                        return false;
                    }
                }
            }
        }
    } else {
        const real_t *r     = xtest[0];
        const real_t *pperp = xtest[1];
        const real_t *ppar  = xtest[2];

        for (len_t i = 0; i < nx1; i++) {
            for (len_t j = 0; j < nx2; j++) {
                for (len_t k = 0; k < nx3; k++) {
                    len_t idx = ((i*nx2 + j)*nx3 + k);
                    real_t Delta = abs(func(r[i], pperp[j], ppar[k]) / dataPP[idx] - 1);

                    if (Delta > TOL) {
                        this->PrintError(
                            "3D interpolation failed for ppar/pperp to ppar/pperp at element ("
                            LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                            "Delta = %e",
                            i, j, k, Delta
                        );
                        return false;
                    }
                }
            }
        }

        const real_t *xi = xtest[1];
        const real_t *p  = xtest[2];

        for (len_t i = 0; i < nx1; i++) {
            for (len_t j = 0; j < nx2; j++) {
                for (len_t k = 0; k < nx3; k++) {
                    len_t idx = ((i*nx2 + j)*nx3 + k);
                    /*real_t p  = sqrt(ppar[k]*ppar[k] + pperp[j]*pperp[j]);
                    real_t xi = ppar[i] / p;

                    real_t Delta = abs(func(r[i], xi, p) / dataPXI[idx] - 1);*/
                    real_t ppar  = p[k]*xi[j];
                    real_t pperp = sqrt(p[k]*p[k] - ppar*ppar);

                    real_t Delta = abs(func(r[i], pperp, ppar) / dataPXI[idx] - 1);

                    if (Delta > TOL) {
                        this->PrintError(
                            "3D interpolation failed for ppar/pperp to p/xi at element ("
                            LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT ", " LEN_T_PRINTF_FMT "). "
                            "Delta = %e",
                            i, j, k, Delta
                        );
                        return false;
                    }
                }
            }
        }
    }

    delete [] dataPP;
    delete [] dataPXI;

    return true;
}

/**
 * Test interpolation method for a general given
 * function, using the specified interpolation method.
 * Note that this function expects that selected interpolation
 * method can reproduce the given function to machine precision.
 *
 * meth: Interpolation method to use.
 * func: Test function to evaluate.
 */
bool Interpolator3D::TestInterpolation3D(
    enum DREAM::FVM::Interpolator3D::interp_method meth,
    function<real_t(real_t, real_t, real_t)> func
) {
    bool success = true;

    struct gridlimits limitsPXI = {0, 1, -1, 1, 0, 1};
    struct gridlimits limitsPP  = {0, 1, 0, 1, -1, 1};

    // Generate (x1,x2,x3) grid and data
    const len_t nx1 = 10, nx2 = 11, nx3 = 12;
    struct griddata *gdPXI = GenerateData(limitsPXI, nx1, nx2, nx3, func);

    DREAM::FVM::Interpolator3D *interpPXI = new DREAM::FVM::Interpolator3D(
        gdPXI->nx1, gdPXI->nx2, gdPXI->nx3,
        gdPXI->x1, gdPXI->x2, gdPXI->x3, gdPXI->y,
        DREAM::FVM::Interpolator3D::GRID_PXI, meth
    );

    success &= EvalInterpolator3D(interpPXI, func);

    delete interpPXI;
    delete gdPXI;

    struct griddata *gdPP = GenerateData(limitsPP, nx1, nx2, nx3, func);

    DREAM::FVM::Interpolator3D *interpPP = new DREAM::FVM::Interpolator3D(
        gdPP->nx1, gdPP->nx2, gdPP->nx3,
        gdPP->x1, gdPP->x2, gdPP->x3, gdPP->y,
        DREAM::FVM::Interpolator3D::GRID_PPARPPERP, meth
    );

    success &= EvalInterpolator3D(interpPP, func);

    delete interpPP;
    delete gdPP;

    return success;
}

/**
 * Test 3D interpolation using the 'nearest' method.
 */
/*bool Interpolator3D::TestNearest() {
    return TestInterpolation3D(
        DREAM::FVM::Interpolator3D::INTERP_NEAREST,
        [](real_t x1, real_t x2, real_t x3) {
            return (4.0*x1 + 2.1*x2 - 7.3*x3 + 2.0);
        }
    );
}*/

/**
 * Test 3D interpolation using the 'linear' method.
 */
bool Interpolator3D::TestLinear() {
    return TestInterpolation3D(
        DREAM::FVM::Interpolator3D::INTERP_LINEAR,
        [](real_t x1, real_t x2, real_t x3) {
            return (4.0*x1 + 2.1*x2 - 7.3*x3 + 2.0);
        }
    );
}

/**
 * Run this test.
 */
bool Interpolator3D::Run(bool) {
    bool success = true;

    /*if (TestNearest())
        this->PrintOK("3D interpolation with 'nearest' method works.");
    else {
        this->PrintError("3D interpolation with 'nearest' method failed.");
        success = false;
    }*/

    if (TestLinear()) {
        this->PrintOK("3D interpolation with 'linear' method works.");
    } else {
        this->PrintError("3D interpolation with 'linear' method failed.");
        success = false;
    }

    return success;
}

