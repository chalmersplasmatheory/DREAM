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
    a->SetCoefficients(this->fr, this->f1, this->f2, this->f1pSqAtZero);
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
 * Returns the delta coefficient given a Peclet number Pe,
 * going smoothly between a central difference scheme in diffusion
 * dominated cases to first-order upwind difference in advection
 * dominated cases (the later characterised by |Pe|>>1) 
 * Solutions are quite sensitive to Pe_threshold, and
 * there is no recipe for how to choose it in general. 
 * Somewhere between 2 and 100 seems good in some cases.
 */
real_t UpwindDelta(real_t A, real_t D, real_t dx){
    real_t Pe_threshold = 10;
    real_t delta_central = 0.5;        
    real_t delta_upwind = (A<0)  - (A>0);

    if(D==0) // upwind for pure advection
        return delta_upwind;

    real_t Pe = dx*A/D;

    real_t limiter;
    real_t x = abs(Pe)/Pe_threshold;
    if(x<=1)
        limiter = 1;
    else
        limiter = 1 / (1 + (x-1)*(x-1)/(x*sqrt(x)) );
    return limiter*delta_central + (1-limiter)*delta_upwind;
}

/** 
 * Set interpolation coefficients dynamically based on the Peclet number
 * of the advection and diffusion coefficients currently stored in memory.
 * Upwind difference for |Pe|>Pe_threshold and central otherwise.
 * Pe_threshold should be at least two, below which central difference 
 * is known to always be stable.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientValuesUpwind() {
    real_t fr, f1, f2, drr, d11, d22;
    real_t DeltaR, DeltaP1, DeltaP2;
    for (len_t ir = 0; ir < this->AdvectionTerm::nr; ir++) {
        const len_t n2 = this->AdvectionTerm::n2[ir];
        const len_t n1 = this->AdvectionTerm::n1[ir];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                fr  = this->AdvectionTerm::fr[ir][j*n1+i] ;
                f1  = this->AdvectionTerm::f1[ir][j*(n1+1)+i] ;
                f2  = this->AdvectionTerm::f2[ir][j*n1+i] ;
                drr = this->DiffusionTerm::drr[ir][j*n1+i];
                d11 = this->DiffusionTerm::d11[ir][j*(n1+1)+i];
                d22 = this->DiffusionTerm::d22[ir][j*n1+i];
                DeltaR  = this->AdvectionTerm::grid->GetRadialGrid()->GetDr(ir);
                DeltaP1 = this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp1(i);
                DeltaP2 = this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp2(j);

                this->AdvectionTerm::deltar[ir][j*n1 + i] = UpwindDelta(fr,drr,DeltaR);
                this->AdvectionTerm::delta1[ir][j*n1 + i] = UpwindDelta(f1,d11,DeltaP1);
                this->AdvectionTerm::delta2[ir][j*n1 + i] = UpwindDelta(f2,d22,DeltaP2);
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

