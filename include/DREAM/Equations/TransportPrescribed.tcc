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
    
    // Interpolate input coefficient onto 'grid'...
    InterpolateCoefficient();
}

/**
 * Destructor.
 */
template<typename T>
DREAM::TransportPrescribed<T>::~TransportPrescribed() {}

/**
 * Function called whenever the computational grid has been
 * rebuilt. Here, we will need to interpolate the prescribed
 * diffusion coefficient onto the new grid.
 */
template<typename T>
bool DREAM::TransportPrescribed<T>::GridRebuilt() {
    constexpr bool va = std::is_same_v<T, DREAM::FVM::AdvectionTerm>;
    constexpr bool vd = std::is_same_v<T, DREAM::FVM::DiffusionTerm>;

    if constexpr (va)
        this->DREAM::FVM::AdvectionTerm::GridRebuilt();
    else if constexpr (vd)
        this->DREAM::FVM::DiffusionTerm::GridRebuilt();
    else
        static_assert(!va && !vd, "The 'TransportPrescribed' term can only be used with 'AdvectionTerm' and 'DiffusionTerm'.");

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
    const len_t N  = this->grid->GetNCells();

    newdata[0] = new real_t[nt*N];

    for (len_t i = 0; i < nt; i++) {
        if (i > 0)
            newdata[i] = newdata[i-1] + N;

        DREAM::FVM::Interpolator3D intp3(
            nr, np1, np2, r, p1, p2, coeff[i],
            momtype, interpmethod, false
        );
        intp3.Eval(this->grid, this->gridtype, newdata[i]);
    }

    if (this->prescribedCoeff != nullptr)
        delete this->prescribedCoeff;

    this->prescribedCoeff = new DREAM::FVM::Interpolator1D(
        nt, N, t, newdata[0]
    );
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
    
    for (len_t ir = 0, offset = 0; ir < nr; ir++) {
        const len_t N = this->grid->GetMomentumGrid(ir)->GetNCells();

        for (len_t j = 0; j < N; j++) {
            constexpr bool va = std::is_same_v<T, DREAM::FVM::AdvectionTerm>;
            constexpr bool vd = std::is_same_v<T, DREAM::FVM::DiffusionTerm>;
            // Set advection/diffusion coefficient...
            if constexpr (va)
                this->fr[ir][j] = c[offset + j];
            else if constexpr (vd)
                this->drr[ir][j] = c[offset + j];
            else
                static_assert(!va && !vd, "The 'TransportPrescribed' term can only be used with 'AdvectionTerm' and 'DiffusionTerm'.");
        }

        offset += N;
    }
}

