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
    a->SetInterpolationCoefficients(this->deltar, this->delta1, this->delta2);
    advectionterms.push_back(a);
}

/**
 * Add a diffusion term to this term.
 */
void AdvectionDiffusionTerm::Add(DiffusionTerm *d) {
    d->SetCoefficients(
        this->drr, this->d11, this->d12, this->d21, this->d22
    );
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
 * Reset the advection and diffusion coefficients of this term.
 */
void AdvectionDiffusionTerm::ResetCoefficients() {
    if (this->advectionterms.size() > 0)
        this->AdvectionTerm::ResetCoefficients();
    if (this->diffusionterms.size() > 0)
        this->DiffusionTerm::ResetCoefficients();

}

/**
 * Rebuild this equation term.
 *
 * t: Simulation time to rebuild term for.
 */
void AdvectionDiffusionTerm::Rebuild(const real_t t, const real_t dt, UnknownQuantityHandler *uqty) {
    this->ResetCoefficients();

    // Rebuild advection-diffusion coefficients
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++){
        (*it)->Rebuild(t, dt, uqty);
    }

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++){
        (*it)->Rebuild(t, dt, uqty);
    }

    // Rebuild interpolation coefficients
    RebuildInterpolationCoefficients(uqty);
}

/**
 * Rebuild the interpolation coefficients.
 */
void AdvectionDiffusionTerm::RebuildInterpolationCoefficients(UnknownQuantityHandler* /*unknowns*/) {
    switch (this->interpolationMethod) {
        case AD_INTERP_CENTRED:  SetInterpolationCoefficientValues(0.5); break;
        case AD_INTERP_BACKWARD: SetInterpolationCoefficientValues(0); break;
        case AD_INTERP_FORWARD:  SetInterpolationCoefficientValues(1); break;
        case AD_INTERP_UPWIND:   SetInterpolationCoefficientValuesUpwind(); break;

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
 * Set interpolation coefficients dynamically based on the Peclet number
 * of the advection and diffusion coefficients currently stored in memory.
 * Upwind difference for |Pe|>Pe_threshold and central otherwise.
 * Pe_threshold should be a maximum of 2, at which point central difference 
 * becomes unstable.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientValuesUpwind() {

    real_t Pe_threshold = 1.8;
    real_t Pe;
    real_t fr, f1, f2, drr, d11, d22;
    real_t d_r, d_1, d_2;
    for (len_t ir = 0; ir < this->AdvectionTerm::nr; ir++) {
        const len_t n2 = this->AdvectionTerm::n2[ir];
        const len_t n1 = this->AdvectionTerm::n1[ir];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                fr = this->AdvectionTerm::fr[ir][j*n1+i] ;
                f1 = this->AdvectionTerm::f1[ir][j*n1+i] ;
                f2 = this->AdvectionTerm::f2[ir][j*n1+i] ;
                drr = this->DiffusionTerm::drr[ir][j*n1+i];
                d11 = this->DiffusionTerm::d11[ir][j*n1+i];
                d22 = this->DiffusionTerm::d22[ir][j*n1+i];
                if(drr==0)
                    d_r = (fr<0) - (fr>0);
                else{
                    Pe = (fr/drr)*this->AdvectionTerm::grid->GetRadialGrid()->GetDr(ir);
                    d_r = 0.5*(1 - (Pe>Pe_threshold) + (Pe<Pe_threshold) );
                }
                if(d11==0)
                    d_1 = (f1<0) - (f1>0);
                else{
                    Pe = (f1/d11)*this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp1(i);
                    d_1 = 0.5*(1 - (Pe>Pe_threshold) + (Pe<Pe_threshold) );
                }
                if(d22==0)
                    d_2 = (f2<0) - (f2>0);
                else{
                    Pe = (f2/d22)*this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp2(i);
                    d_2 = 0.5*(1 - (Pe>Pe_threshold) + (Pe<Pe_threshold) );
                }                
                this->AdvectionTerm::deltar[ir][j*n1 + i] = d_r;
                this->AdvectionTerm::delta1[ir][j*n1 + i] = d_1;
                this->AdvectionTerm::delta2[ir][j*n1 + i] = d_2;
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

