/**
 * Implementation of a 3D interpolation class.
 */

#include <cmath>
#include <gsl/gsl_interp.h>
#include "FVM/Interpolator3D.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 */
Interpolator3D::Interpolator3D(
    const len_t nx1, const len_t nx2, const len_t nx3,
    const real_t *x1, const real_t *x2, const real_t *x3,
    const real_t *y,
    enum momentumgrid_type type, enum interp_method meth
) : nx1(nx1), nx2(nx2), nx3(nx3),
    x1(x1), x2(x2), x3(x3), y(y),
    gridtype(type), method(meth) {

    this->acc1 = gsl_interp_accel_alloc();
    this->acc2 = gsl_interp_accel_alloc();
    this->acc3 = gsl_interp_accel_alloc();
}

/**
 * Destructor.
 */
Interpolator3D::~Interpolator3D() {
    gsl_interp_accel_free(this->acc3);
    gsl_interp_accel_free(this->acc2);
    gsl_interp_accel_free(this->acc1);

    delete [] this->y;
    delete [] this->x3;
    delete [] this->x2;
    delete [] this->x1;
}

/**
 * Evaluate the data of this object on the
 * given computational grid.
 *
 * grid: Computational grid to evaluate the interpolator object on.
 * type: Type of the momentum part of 'grid' (i.e. p/xi or ppar/pperp).
 */
const real_t *Interpolator3D::Eval(FVM::Grid *grid, enum momentumgrid_type type) {
    // XXX Here we assume that all momentum grids are the same
    return this->Eval(
        grid->GetNr(),
        grid->GetMomentumGrid(0)->GetNp1(),
        grid->GetMomentumGrid(0)->GetNp2(),
        grid->GetRadialGrid()->GetR(),
        grid->GetMomentumGrid(0)->GetP1(),
        grid->GetMomentumGrid(0)->GetP2(),
        type
    );
}

/**
 * Evaluate the data of this object on the
 * given computational grid.
 *
 * nx1, nx2, nx3: Length of each dimension of the grid to evaluate
 *                the data on.
 * x1, x2, x3:    Coordinates of the grid to evaluate the data on.
 * type:          Type of the momentum part of the grid (i.e. p/xi
 *                or ppar/pperp).
 */
const real_t *Interpolator3D::Eval(
    const len_t nx1, const len_t nx2, const len_t nx3,
    const real_t *x1, const real_t *x2, const real_t *x3,
    enum momentumgrid_type type
) {
    real_t *data = new real_t[nx1*nx2*nx3];
    const enum interp_method meth = this->method;

    #define EVAL(X1,X2,X3) \
        for (len_t k = 0; k < nx1; k++) { \
            for (len_t j = 0; j < nx2; j++) { \
                for (len_t i = 0; i < nx3; i++) { \
                    const len_t idx = (k*nx2 + j)*nx3 + i; \
                    if (meth == INTERP_NEAREST) { \
                        data[idx] = this->_eval_nearest((X1), (X2), (X3)); \
                    } else { \
                        data[idx] = this->_eval_linear((X1), (X2), (X3)); \
                    } \
                } \
            } \
        }

    if (type == this->gridtype) {
        EVAL(x1[k], x2[j], x3[i]);
    } else if (type == GRID_PXI) {
        const real_t *p  = x3;
        const real_t *xi = x2;
        #define PPAR (p[i]*xi[j])
        #define PPERP (p[i]*sqrt(1-xi[j]*xi[j]))
        EVAL(x1[k], PPERP, PPAR);
        #undef PPERP
        #undef PPAR
    } else if (type == GRID_PPARPPERP) {
        const real_t *ppar  = x3;
        const real_t *pperp = x2;
        #define P (sqrt(ppar[i]*ppar[i] + pperp[j]*pperp[j]))
        #define XI (ppar[i]/P)
        EVAL(x1[k], XI, P);
        #undef XI
        #undef P
    }

    #undef EVAL
    
    return data;
}

/**
 * Evaluate a single point on the grid using the
 * 'nearest' interpolation algorithm.
 */
real_t Interpolator3D::_eval_nearest(
    const real_t x1, const real_t x2, const real_t x3
) {
    len_t ix1 = _find_x1(x1);
    len_t ix2 = _find_x2(x2);
    len_t ix3 = _find_x3(x3);

    #define CORRECT(arr) \
        if ((i ## arr) +1 < (this->n ## arr)) {\
            if (abs(this-> arr [(i ## arr)] - arr) > abs(this-> arr [(i ## arr)+1] - arr)) \
                (i ## arr)++; \
        }

    CORRECT(x1);
    CORRECT(x2);
    CORRECT(x3);

    #undef CORRECT

    return this->y[(ix1*nx2 + ix2)*nx3 + ix3];
}

/**
 * Evaluate a single point on the grid using the
 * 'linear' interpolation algorithm.
 */
real_t Interpolator3D::_eval_linear(
    const real_t x1, const real_t x2, const real_t x3
) {
    len_t ix10 = _find_x1(x1);
    len_t ix20 = _find_x2(x2);
    len_t ix30 = _find_x3(x3);

    if (ix10+1 == this->nx1) ix10--;
    if (ix20+1 == this->nx2) ix20--;
    if (ix30+1 == this->nx3) ix30--;

    len_t ix11 = ix10 + 1;
    len_t ix21 = ix20 + 1;
    len_t ix31 = ix30 + 1;

    // Check for single grid points
    if (ix11 == this->nx1) ix11 = ix10;
    if (ix21 == this->nx2) ix21 = ix20;
    if (ix31 == this->nx3) ix31 = ix30;

    #define IDX(X1,X2,X3) (((X1)*nx2 + (X2))*nx3 + (X3))

    real_t y000 = this->y[IDX(ix10, ix20, ix30)];
    real_t y100 = this->y[IDX(ix11, ix20, ix30)];
    real_t y010 = this->y[IDX(ix10, ix21, ix30)];
    real_t y001 = this->y[IDX(ix10, ix20, ix31)];
    real_t y110 = this->y[IDX(ix11, ix21, ix30)];
    real_t y101 = this->y[IDX(ix11, ix20, ix31)];
    real_t y011 = this->y[IDX(ix10, ix21, ix31)];
    real_t y111 = this->y[IDX(ix11, ix21, ix31)];

    real_t x1d=0, x2d=0, x3d=0;
    if (ix10 != ix11) x1d = (x1-this->x1[ix10]) / (this->x1[ix11] - this->x1[ix10]);
    if (ix20 != ix21) x2d = (x2-this->x2[ix20]) / (this->x2[ix21] - this->x2[ix20]);
    if (ix30 != ix31) x3d = (x3-this->x3[ix30]) / (this->x3[ix31] - this->x3[ix30]);

    real_t y00 = y000*(1 - x1d) + y100*x1d;
    real_t y01 = y001*(1 - x1d) + y101*x1d;
    real_t y10 = y010*(1 - x1d) + y110*x1d;
    real_t y11 = y011*(1 - x1d) + y111*x1d;

    real_t y0 = y00*(1 - x2d) + y10*x2d;
    real_t y1 = y01*(1 - x2d) + y11*x2d;

    return (y0*(1 - x3d) + y1*x3d);
}

/**
 * Locate the index i that is such that
 *
 *   xarr[i] <= x < xarr[i+1]
 *
 * x:    Value of 'x' to look up.
 * nx:   Number of elements in 'xarr'.
 * xarr: Array to search in.
 * acc:  GSL accelerator object used to speed up search.
 */
len_t Interpolator3D::_find_x(
    const real_t x, const len_t nx,
    const real_t *xarr, gsl_interp_accel *acc
) {
    return (len_t)gsl_interp_accel_find(acc, xarr, nx, x);
}

