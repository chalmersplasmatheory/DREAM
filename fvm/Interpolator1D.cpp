/**
 * Implementation of 1D vector interpolation.
 *
 * This object takes two vectors as input -- one vector 'x' with 'nx' elements,
 * and one vector 'y' with ny = nx*nblocks elements -- where y = y(x). The
 * vector y is assumed to consist of 'nblocks' independent blocks of data which
 * all depend on the parameter 'x'. Upon calling 'Eval()', the value of 'y' in
 * each block will be evaluated and returned as a new vector with 'nblocks'
 * elements.
 *
 * We assume that 'x' is sorted. It may however be either monotonously
 * increasing or decreasing.
 *
 * The structure of 'y' should be
 *
 *   y(x0,b0) y(x0,b1) y(x0,b2) ... y(x0,bN) y(x1,b0) ... y(tM,bN)
 *
 * i.e. data belonging to the same time are located close together.
 */

#include <cmath>
#include "FVM/Interpolator1D.hpp"
#include <gsl/gsl_machine.h>


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * nx:      Number of x points.
 * nblocks: Number of blocks in y data.
 * x:       List of parameter values.
 * y:       List of elements to interpolate in.
 */
Interpolator1D::Interpolator1D(
    const len_t nx, const len_t nblocks, const real_t *x, const real_t *y,
    enum interp_method meth, bool owns_data
) : nx(nx), nblocks(nblocks), x(x), y(y), method(meth), owns_data(owns_data) {
    // Since the 'nearest' interpolation method returns an
    // exact copy of some of the data in 'y', we won't need
    // a buffer for that method.
    if (meth == INTERP_LOGARITHMIC){
        this->logy = new real_t[nx*nblocks];
        len_t i, ix, ib;
        for (ix = 1; ix < nx; ix++)
            for (ib = 0; ib < nblocks; ib++){
                i = ix*nblocks + ib;
                if (y[i] > GSL_DBL_MIN)
                    logy[i] = log(y[i]);
                else
                    logy[i] = GSL_LOG_DBL_MIN;
            }
    }
    if (meth != INTERP_NEAREST)
        this->buffer = new real_t[nblocks];

    if (nx > 1)
        this->xIncreasing = (x[1]>x[0]);

    // Verify that x is monotonically increasing/decreasing...
    for (len_t i = 1; i < nx; i++)
        if ((xIncreasing && x[i] <= x[i-1]) ||
            (!xIncreasing && x[i] >= x[i-1]))
                throw Interpolator1DException(
                    "interpolator1d: x data must be monotonically increasing/decreasing."
                );
}

/**
 * Destructor.
 */
Interpolator1D::~Interpolator1D() {
    if (buffer != nullptr)
        delete [] this->buffer;
    if (logy != nullptr)
        delete [] this->logy;
	if (this->owns_data) {
		delete [] this->x;
		delete [] this->y;
	}
}

/**
 * Evaluate 'y' in the given 'x' point.
 *
 * NOTE: The returned pointer is owned by this pointer and
 * should NOT be free'd outside of this object!
 */
const real_t *Interpolator1D::Eval(const real_t x) {
    switch (this->method) {
        case INTERP_NEAREST: return _eval_nearest(x);
        case INTERP_LOGARITHMIC: return _eval_logarithmic(x);
        case INTERP_LINEAR: return _eval_linear(x);

        // We shouldn't end up here, so just return
        // 'nullptr' to make the compiler shut up
        default:
            return nullptr;
    }
}

/**
 * Locate the index of an 'x' that is before the given
 * value of 'x' in the x vector.
 */
len_t Interpolator1D::_find_x(const real_t xv) {
    int a = 0;
    int b = nx-1;

    while (abs(a-b) > 1) {
        int c = (a+b)/2;
        if (x[c] >= xv) {
            if (xIncreasing) b = c;
            else a = c;
        } else {
            if (xIncreasing) a = c;
            else b = c;
        }
    }

    return a;
}

/**
 * Linear interpolation.
 */
const real_t *Interpolator1D::_eval_linear(const real_t xv) {
    len_t ix = _find_x(xv);

    if (ix+1 >= nx) {
        for (len_t i = 0; i < nblocks; i++)
            buffer[i] = y[(nx-1)*nblocks + i];

        return buffer;
    }

    const real_t x1  = x[ix];
    const real_t x2  = x[ix+1];
    const real_t ddx = (xv-x1) / (x2-x1);
    for (len_t i = 0; i < nblocks; i++) {
        const real_t y1 = y[ix*nblocks + i];
        const real_t y2 = y[(ix+1)*nblocks + i];

        buffer[i] = (1-ddx)*y1 + ddx*y2;
    }

    return buffer;
}

/**
 * Logarithmic interpolation.
 */
const real_t *Interpolator1D::_eval_logarithmic(const real_t xv) {
    len_t ix = _find_x(xv);

    if (ix+1 >= nx) {
        for (len_t i = 0; i < nblocks; i++)
            buffer[i] = logy[(nx-1)*nblocks + i];

        return buffer;
    }

    const real_t x1  = x[ix];
    const real_t x2  = x[ix+1];
    const real_t ddx = (xv-x1) / (x2-x1);
    for (len_t i = 0; i < nblocks; i++) {
        const real_t y1 = logy[ix*nblocks + i];
        const real_t y2 = logy[(ix+1)*nblocks + i];

        buffer[i] = exp((1-ddx)*y1 + ddx*y2);
    }

    return buffer;
}

/**
 * Fetch 'y' at the parameter value closest to 'xv'.
 */
const real_t *Interpolator1D::_eval_nearest(const real_t xv) {
    len_t ix = _find_x(xv);

    if (ix+1 >= nx)
        return y+((nx-1)*nblocks);
    else if (fabs(xv-x[ix]) < fabs(xv-x[ix+1]))
        return y+(ix*nblocks);
    else
        return y+((ix+1)*nblocks);
}

