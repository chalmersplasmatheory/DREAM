/**
 * Implementation of a combined advection-diffusion term.
 */

#include <algorithm>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"


using namespace DREAM::FVM;
using namespace std;

/**
 * Add an advection term to this term.
 */
void AdvectionDiffusionTerm::Add(AdvectionTerm *a) {
    a->SetCoefficients(this->fr, this->f1, this->f2);
    advectionterms.push_back(a);
}

/**
 * Add a diffusion term to this term.
 */
void AdvectionDiffusionTerm::Add(DiffusionTerm *d) {
    d->SetCoefficients(this->drr, this->d11, this->d12, this->d21, this->d22);
    diffusionterms.push_back(d);
}

/**
 * Returns the number of non-zero elements per row
 * inserted by this term into a linear operator matrix.
 */
len_t AdvectionDiffusionTerm::GetNumberOfNonZerosPerRow() const {
    len_t nnz = 0;
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    return nnz;
}

/**
 * Returns the number of non-zero elements per row
 * inserted by this term into a jacobian matrix.
 */
len_t AdvectionDiffusionTerm::GetNumberOfNonZerosPerRow_jac() const {
    len_t nnz = 0;
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow_jac());

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow_jac());

    return nnz;
}

/**
 * Rebuild this equation term.
 *
 * t: Simulation time to rebuild term for.
 */
void AdvectionDiffusionTerm::Rebuild(const real_t t, const real_t dt, UnknownQuantityHandler *uqty) {
    // Rebuild advection-diffusion coefficients
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Rebuild interpolation coefficients
    RebuildInterpolationCoefficients();
}

/**
 * Rebuild the interpolation coefficients.
 */
void AdvectionDiffusionTerm::RebuildInterpolationCoefficients() {
    switch (this->interpolationMethod) {
        case AD_INTERP_CENTRED:  SetInterpolationCoefficientValues(0.5); break;
        case AD_INTERP_BACKWARD: SetInterpolationCoefficientValues(0); break;
        case AD_INTERP_FORWARD:  SetInterpolationCoefficientValues(1); break;

        default:
            throw EquationTermException(
                "Unrecognized advection/diffusion interpolation method specified: %d.",
                this->interpolationMethod
            );
    }
}

/**
 * Set all interpolation coefficients to the same value.
 *
 * v: Value to assign to all interpolation coefficients.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientValues(const real_t v) {
    for (len_t k = 0; k < this->AdvectionTerm::nr; k++) {
        const len_t n2 = this->AdvectionTerm::n2[k];
        const len_t n1 = this->AdvectionTerm::n1[k];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                this->AdvectionTerm::deltar[k][j*n1 + i] = v;
                this->AdvectionTerm::delta1[k][j*n1 + i] = v;
                this->AdvectionTerm::delta2[k][j*n1 + i] = v;
            }
        }
    }
}

/**
 * Save the advection and diffusion coefficients of this
 * object to the specified file.
 */
void AdvectionDiffusionTerm::SaveCoefficientsSFile(const string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    this->SaveCoefficientsSFile(sf);
    sf->Close();
}
void AdvectionDiffusionTerm::SaveCoefficientsSFile(SFile *sf) {
    this->AdvectionTerm::SaveCoefficientsSFile(sf);
    this->DiffusionTerm::SaveCoefficientsSFile(sf);
}

