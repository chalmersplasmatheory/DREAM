/**
 * Implementation of a diffusive transport term which can be applied to both
 * kinetic and fluid grids, and which allows one to prescribe the diffusion
 * coefficient in time and phase space.
 */

#include <type_traits>
#include "DREAM/Equations/TransportPrescribed.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Interpolator3D.hpp"


/**
 * Constructor.
 */
template<typename T>
DREAM::TransportPrescribed<T>::TransportPrescribed(
    DREAM::FVM::Grid *grid,
    const len_t nt, const len_t nr, const len_t np1, const len_t np2,
    const real_t **coeff, const real_t *t, const real_t *r,
    const real_t *p1, const real_t *p2, enum DREAM::FVM::Interpolator3D::momentumgrid_type inptype,
    enum DREAM::FVM::Interpolator3D::momentumgrid_type gridtype,
    enum DREAM::FVM::Interpolator3D::interp_method interpmethod,
    bool allocCoefficients
) : T(grid, allocCoefficients),
    nt(nt), nr(nr), np1(np1), np2(np2),
    coeff(coeff), t(t), r(r), p1(p1), p2(p2),
    momtype(inptype), gridtype(gridtype), interpmethod(interpmethod) {

    this->T::SetName("TransportPrescribed");
    
    // Interpolate input coefficient onto 'grid'...
    InterpolateCoefficient();
}

/**
 * Destructor.
 */
namespace DREAM {
template<typename T>
TransportPrescribed<T>::~TransportPrescribed() {
    if (this->prescribedCoeff != nullptr)
        delete this->prescribedCoeff;
    if (this->interpolateddata != nullptr) {
        delete [] this->interpolateddata;
    }

    delete [] this->coeff[0];
    delete [] this->coeff;
    delete [] this->t;
    delete [] this->r;
    delete [] this->p1;
    delete [] this->p2;
}
}


/**
 * Function called whenever the computational grid has been
 * rebuilt. Here, we will need to interpolate the prescribed
 * diffusion coefficient onto the new grid.
 */
template<typename T>
bool DREAM::TransportPrescribed<T>::GridRebuilt() {
    //  Call GridRebuilt() on base class...
    this->T::GridRebuilt();

    // Interpolate prescribed coefficient onto new
    // phase space grid...
    InterpolateCoefficient();

    return true;
}

/**
 * Interpolate the input coefficient onto the phase space
 * grid used by this EquationTerm.
 */
template<typename T>
void DREAM::TransportPrescribed<T>::InterpolateCoefficient() {
    real_t **newdata = new real_t*[nt];
    // Drr is defined on the radial flux grid so we
    // XXX assume that all momentum grids are the same
    // and one momentum grid's worth of cells to the total
    // number of cells...
    const len_t N  =
        this->grid->GetNCells() + this->grid->GetMomentumGrid(0)->GetNCells();

    newdata[0] = new real_t[nt*N];
    real_t* newt = new real_t[nt];

    for (len_t i = 0; i < nt; i++) {
        if (i > 0)
            newdata[i] = newdata[i-1] + N;
        newt[i] = t[i];

        DREAM::FVM::Interpolator3D intp3(
            nr, np2, np1, r, p2, p1, coeff[i],
            momtype, interpmethod, false
        );
        intp3.Eval(this->grid, this->gridtype, FVM::FLUXGRIDTYPE_RADIAL, newdata[i]);
    }

    if (this->prescribedCoeff != nullptr) {
        delete this->prescribedCoeff;
        delete [] this->interpolateddata;
    }

    this->prescribedCoeff = new DREAM::FVM::Interpolator1D(
        nt, N, newt, newdata[0]
    );

    // This data is now used by 'prescribedCoeff', but we
    // need to keep a pointer in this class so that we can
    // clean it up later...
    this->interpolateddata = newdata;
}

/**
 * Rebuild this term by evaluating and setting the diffusion
 * coefficient for the next time step.
 */
template<typename T>
void DREAM::TransportPrescribed<T>::Rebuild(
    const real_t t, const real_t, DREAM::FVM::UnknownQuantityHandler*
) {
    const real_t *c = this->prescribedCoeff->Eval(t);
    const len_t nr = this->grid->GetNr();
    // XXX here we assume that all momentum grids are the same...
    const len_t N = this->grid->GetMomentumGrid(0)->GetNCells();
    
    // Iterate over the radial flux grid...
    for (len_t ir = 0, offset = 0; ir < nr+1; ir++) {
        for (len_t j = 0; j < N; j++) {
            this->_setcoeff(ir, j, c[offset+j]);
        }

        offset += N;
    }
}

