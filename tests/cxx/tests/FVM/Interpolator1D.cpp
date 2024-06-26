/**
 * Test for the 1D interpolation object.
 */

#include "FVM/Interpolator1D.hpp"
#include "Interpolator1D.hpp"


using namespace DREAMTESTS::FVM;
using namespace std;


/**
 * Implementation of an interpolation test with a linear
 * function as test function. The input function 'cmp' does
 * the following:
 *
 * Check if the given interpolated array 'intp' was
 * calculated correctly, by comparing to the exact function
 * 'f'. The input array consists of 'nblocks' blocks.
 * 
 * nblocks: Number of elements in 'intp'
 * intp:    Array containing interpolated values
 * a:       Input parameter to 'f'
 * f:       Function to compare against
 */
bool Interpolator1D::TestLinear_general(
    enum DREAM::FVM::Interpolator1D::interp_method meth,
    function<bool(
        len_t, len_t, real_t, const real_t*, const real_t*, const real_t*,
        const function<real_t(const real_t, const real_t)>&)> cmp
) {
    bool success = true;
    const len_t nx = 10, nblocks = 13;
    const len_t ntests = 53;

    // Test input
    const real_t xtest[ntests] = {
        0.1715,0.0677,0.4803,0.9990,0.8157,0.1485,0.7714,0.9366,0.6683,0.5587,
        0.0241,0.9430,0.0142,0.7211,0.4934,0.9856,0.2388,0.2387,0.9087,0.1552,
        0.3993,0.9062,0.4599,0.2206,0.9581,0.8101,0.6662,0.0770,0.6825,0.2404,
        0.0278,0.1359,0.2677,0.7440,0.2253,0.1133,0.1772,0.7712,0.6218,0.7722,
        0.1514,0.1462,0.3076,0.0934,0.1439,0.2037,0.2709,0.9787,0.9291,0.5303,
        0.8650,0.4277,0.4109
    };

    function<real_t(const real_t, const real_t)> f =
        [](const real_t x, const real_t a) { return a*x + 5; };

    // Generate input parameters
    real_t *x = new real_t[nx];
    real_t *a = new real_t[nblocks];
    for (len_t i = 0; i < nx; i++)
        x[i] = ((real_t)i) / (nx-1);
    for (len_t i = 0; i < nblocks; i++)
        a[i] = -1.5 + 4*((real_t)i) / (nblocks-1);

    // Generate data
    real_t *v = new real_t[nx*nblocks];
    for (len_t j = 0; j < nx; j++)
        for (len_t i = 0; i < nblocks; i++)
            v[j*(nblocks) + i] = f(x[j], a[i]);

    auto *intp = new DREAM::FVM::Interpolator1D(nx, nblocks, x, v, meth, false);
    for (len_t i = 0; i < ntests && success; i++) {
        const real_t *b = intp->Eval(xtest[i]);

        if (!cmp(nblocks, nx, xtest[i], x, b, a, f)) {
            this->PrintError(
                "Linear 1D interpolation failed at test value index " LEN_T_PRINTF_FMT
                ".", i
            );
            success = false;
        }
    }

    delete intp;

    // Test with reversed x array
    for (len_t i = 0; i < nx; i++)
        x[i] = 1-((real_t)i) / (nx-1);

    for (len_t j = 0; j < nx; j++)
        for (len_t i = 0; i < nblocks; i++)
            v[j*(nblocks) + i] = f(x[j], a[i]);

    intp = new DREAM::FVM::Interpolator1D(nx, nblocks, x, v, meth, false);
    for (len_t i = 0; i < ntests && success; i++) {
        const real_t *b = intp->Eval(xtest[i]);

        if (!cmp(nblocks, nx, xtest[i], x, b, a, f)) {
            this->PrintError(
                "Linear 1D interpolation failed at test value index " LEN_T_PRINTF_FMT
                ".", i
            );
            success = false;
        }
    }

    delete intp;
    delete [] v;
    delete [] a;
    delete [] x;

    return success;
}

/**
 * Test 'nearest' interpolation.
 */
bool Interpolator1D::TestLinear() {
    return TestLinear_general(
        DREAM::FVM::Interpolator1D::INTERP_LINEAR, [](
        len_t nblocks, len_t /*nx*/, real_t x, const real_t*, const real_t *intp,
        const real_t *a, const function<real_t(const real_t, const real_t)>& f
    ) {
        const real_t TOL = 10 * std::numeric_limits<real_t>::epsilon();

        for (len_t i = 0; i < nblocks; i++) {
            real_t fval = f(x,a[i]);
            if (abs(intp[i]/fval - 1) > TOL)
                return false;
        }

        return true;
    });
}

/**
 * Test 'nearest' interpolation.
 */
bool Interpolator1D::TestNearest() {
    return TestLinear_general(
        DREAM::FVM::Interpolator1D::INTERP_NEAREST, [](
        len_t nblocks, len_t nx, real_t x, const real_t *xvec, const real_t *intp,
        const real_t *a, const function<real_t(const real_t, const real_t)>& f
    ) {
        const real_t TOL = 10*std::numeric_limits<real_t>::epsilon();

        for (len_t i = 0; i < nblocks; i++) {
            // Locate x-value (in a slow way)
            len_t ix = 0;
            while (ix+1 < nx && abs(x-xvec[ix+1]) < abs(x-xvec[ix]))
                ix++;

            real_t fval = f(xvec[ix], a[i]);
            if (abs(intp[i]/fval - 1) > TOL)
                return false;
        }

        return true;
    });
}

/**
 * Run this test.
 */
bool Interpolator1D::Run(bool) {
    bool success = true;
    if ((success=TestNearest()))
        this->PrintOK("1D interpolation with 'nearest' method works.");
    else
        this->PrintError("1D interpolation with 'nearest' method failed.");

    if ((success=TestLinear()))
        this->PrintOK("1D interpolation with 'linear' method works.");
    else
        this->PrintError("1D interpolation with 'linear' method failed.");

    return success;
}

