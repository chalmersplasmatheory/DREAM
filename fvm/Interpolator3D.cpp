/**
 * Implementation of a 3D interpolation class.
 */

#include <cmath>
#include <gsl/gsl_interp.h>
#include "FVM/FVMException.hpp"
#include "FVM/Interpolator3D.hpp"
#include <gsl/gsl_machine.h>

using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 */
Interpolator3D::Interpolator3D(
    const len_t nx1, const len_t nx2, const len_t nx3,
    const real_t *x1, const real_t *x2, const real_t *x3,
    const real_t *y,
    enum momentumgrid_type type, enum interp_method meth,
    bool ownsArrays
) : nx1(nx1), nx2(nx2), nx3(nx3),
    x1(x1), x2(x2), x3(x3), y(y),
    gridtype(type), method(meth),
    ownsArrays(ownsArrays) {
    
    if (meth == INTERP_LOGARITHMIC){
        this->logy = new real_t[nx1*nx2*nx3];
        real_t EXPLOGMIN = exp(GSL_LOG_DBL_MIN);
        len_t i, i1, i2, i3;
        for (i1 = 0; i1 < nx1; i1++)
            for (i2 = 0; i2 < nx2; i2++)
                for (i3 = 0; i3 < nx3; i3++){
                    i = ((i1)*nx2 + (i2))*nx3 + (i3);
                    if (y[i] > EXPLOGMIN)
                        logy[i] = log(y[i]);
                    else
                        logy[i] = GSL_LOG_DBL_MIN;
                }
    }

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


    if (logy != nullptr)
        delete [] this->logy;

    if (this->ownsArrays) {
        delete [] this->y;

        if (this->x3 != nullptr)
            delete [] this->x3;
        if (this->x2 != nullptr)
            delete [] this->x2;
        if (this->x1 != nullptr)
            delete [] this->x1;
    }
}

/**
 * Evaluate the data of this object on the
 * given computational grid.
 *
 * grid: Computational grid to evaluate the interpolator object on.
 * type: Type of the momentum part of 'grid' (i.e. p/xi or ppar/pperp).
 * fgt:  Type of grid (either the distribution grid, r flux grid, p1
 *       flux grid or p2 flux grid).
 * out:  Array to store interpolated data in. If 'nullptr', new memory
 *       is allocated and must later be deleted by the caller.
 */
const real_t *Interpolator3D::Eval(
    FVM::Grid *grid, enum momentumgrid_type type,
    enum fluxGridType fgt, real_t *out
) {
    // XXX Here we assume that all momentum grids are the same
    switch (fgt) {
        case FLUXGRIDTYPE_DISTRIBUTION:
            return this->Eval(
                grid->GetNr(),
                grid->GetMomentumGrid(0)->GetNp2(),
                grid->GetMomentumGrid(0)->GetNp1(),
                grid->GetRadialGrid()->GetR(),
                grid->GetMomentumGrid(0)->GetP2(),
                grid->GetMomentumGrid(0)->GetP1(),
                type, out
            );

        case FLUXGRIDTYPE_RADIAL:
            return this->Eval(
                grid->GetNr()+1,
                grid->GetMomentumGrid(0)->GetNp2(),
                grid->GetMomentumGrid(0)->GetNp1(),
                grid->GetRadialGrid()->GetR_f(),
                grid->GetMomentumGrid(0)->GetP2(),
                grid->GetMomentumGrid(0)->GetP1(),
                type, out
            );

        case FLUXGRIDTYPE_P2:
            return this->Eval(
                grid->GetNr(),
                grid->GetMomentumGrid(0)->GetNp2()+1,
                grid->GetMomentumGrid(0)->GetNp1(),
                grid->GetRadialGrid()->GetR(),
                grid->GetMomentumGrid(0)->GetP2_f(),
                grid->GetMomentumGrid(0)->GetP1(),
                type, out
            );

        case FLUXGRIDTYPE_P1:
            return this->Eval(
                grid->GetNr(),
                grid->GetMomentumGrid(0)->GetNp2(),
                grid->GetMomentumGrid(0)->GetNp1()+1,
                grid->GetRadialGrid()->GetR(),
                grid->GetMomentumGrid(0)->GetP2(),
                grid->GetMomentumGrid(0)->GetP1_f(),
                type, out
            );

        default:
            throw FVMException("Unrecognized flux grid type specified: %d.", fgt);
    }
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
 * out:           Array to store interpolated data in. If 'nullptr',
 *                new memory is allocated and must later be deleted
 *                by the caller.
 */
const real_t *Interpolator3D::Eval(
    const len_t nx1, const len_t nx2, const len_t nx3,
    const real_t *x1, const real_t *x2, const real_t *x3,
    enum momentumgrid_type type, real_t *out
) {
    if (out == nullptr)
        out = new real_t[nx1*nx2*nx3];

    const enum interp_method meth = this->method;

    #define EVAL(X1,X2,X3) \
        for (len_t k = 0; k < nx1; k++) { \
            for (len_t j = 0; j < nx2; j++) { \
                for (len_t i = 0; i < nx3; i++) { \
                    const len_t idx = (k*nx2 + j)*nx3 + i; \
                    if (meth == INTERP_NEAREST) { \
                        out[idx] = this->_eval_nearest((X1), (X2), (X3)); \
                    } else if (meth == INTERP_LOGARITHMIC) { \
                        out[idx] = this->_eval_logarithmic((X1), (X2), (X3)); \
                    } else { \
                        out[idx] = this->_eval_linear((X1), (X2), (X3)); \
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
    
    return out;
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

    #define CORRECT(arr) do { \
        if ((i ## arr) +1 < (this->n ## arr)) {\
            if (abs(this-> arr [(i ## arr)] - arr) > abs(this-> arr [(i ## arr)+1] - arr)) \
                (i ## arr)++; \
        }} while (false)

    if (nx1 > 1) CORRECT(x1);
    if (nx2 > 1) CORRECT(x2);
    if (nx3 > 1) CORRECT(x3);

    #undef CORRECT

    return this->y[(ix1*nx2 + ix2)*nx3 + ix3];
}

/**
 * Evaluate a single point on the grid using the
 * 'logarithmic' interpolation algorithm.
 */
real_t Interpolator3D::_eval_logarithmic(
    const real_t x1, const real_t x2, const real_t x3
) {
    len_t ix10 = _find_x1(x1);
    len_t ix20 = _find_x2(x2);
    len_t ix30 = _find_x3(x3);

    if (this->nx1 > 1 && ix10+1 == this->nx1) ix10--;
    if (this->nx2 > 1 && ix20+1 == this->nx2) ix20--;
    if (this->nx3 > 1 && ix30+1 == this->nx3) ix30--;

    len_t ix11 = ix10 + 1;
    len_t ix21 = ix20 + 1;
    len_t ix31 = ix30 + 1;

    // Check for single grid points
    if (ix11 == this->nx1) ix11 = ix10;
    if (ix21 == this->nx2) ix21 = ix20;
    if (ix31 == this->nx3) ix31 = ix30;

    #define IDX(X1,X2,X3) (((X1)*nx2 + (X2))*nx3 + (X3))

    real_t y000 = this->logy[IDX(ix10, ix20, ix30)];
    real_t y100 = this->logy[IDX(ix11, ix20, ix30)];
    real_t y010 = this->logy[IDX(ix10, ix21, ix30)];
    real_t y001 = this->logy[IDX(ix10, ix20, ix31)];
    real_t y110 = this->logy[IDX(ix11, ix21, ix30)];
    real_t y101 = this->logy[IDX(ix11, ix20, ix31)];
    real_t y011 = this->logy[IDX(ix10, ix21, ix31)];
    real_t y111 = this->logy[IDX(ix11, ix21, ix31)];

    real_t x1d=0, x2d=0, x3d=0;
    if (this->x1 != nullptr)
        if (ix10 != ix11) x1d = (x1-this->x1[ix10]) / (this->x1[ix11] - this->x1[ix10]);
    if (this->x2 != nullptr)
        if (ix20 != ix21) x2d = (x2-this->x2[ix20]) / (this->x2[ix21] - this->x2[ix20]);
    if (this->x3 != nullptr)
        if (ix30 != ix31) x3d = (x3-this->x3[ix30]) / (this->x3[ix31] - this->x3[ix30]);

    real_t y00 = y000*(1 - x1d) + y100*x1d;
    real_t y01 = y001*(1 - x1d) + y101*x1d;
    real_t y10 = y010*(1 - x1d) + y110*x1d;
    real_t y11 = y011*(1 - x1d) + y111*x1d;

    real_t y0 = y00*(1 - x2d) + y10*x2d;
    real_t y1 = y01*(1 - x2d) + y11*x2d;

    return exp(y0*(1 - x3d) + y1*x3d);
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

    if (this->nx1 > 1 && ix10+1 == this->nx1) ix10--;
    if (this->nx2 > 1 && ix20+1 == this->nx2) ix20--;
    if (this->nx3 > 1 && ix30+1 == this->nx3) ix30--;

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
    if (this->x1 != nullptr)
        if (ix10 != ix11) x1d = (x1-this->x1[ix10]) / (this->x1[ix11] - this->x1[ix10]);
    if (this->x2 != nullptr)
        if (ix20 != ix21) x2d = (x2-this->x2[ix20]) / (this->x2[ix21] - this->x2[ix20]);
    if (this->x3 != nullptr)
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
    if (nx == 1) return 0;
    else return (len_t)gsl_interp_accel_find(acc, xarr, nx, x);
}

