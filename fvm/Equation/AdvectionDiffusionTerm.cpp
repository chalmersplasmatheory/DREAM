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


    /**
     * Add additional damping of the flux limiter schemes 
     * as iteration in non-linear solver becomes large.
     * Adds an additional damping factor that linearly goes from 1 
     * at iteration=itThresh to minDamping at iteration=itMax.
     */
    if(withDynamicFluxLimiterDamping){
        real_t minDamping = 0.8;
        len_t itMax = 100;
        len_t itThresh = 30;

        // Manually count iteration number indirectly
        if((this->t_prev==t) && (this->dt_prev==dt))
            iteration +=1;
        else {
            this->iteration = 1;
            dampingWithIteration = 1.0;
        }
        if(this->iteration>itThresh)
            dampingWithIteration = std::max(minDamping, 
                1.0 - ((1.0-minDamping)*(iteration-itThresh))/(itMax-itThresh));
        this->t_prev = t;
        this->dt_prev = dt;
    }

    // Rebuild interpolation coefficients
    RebuildInterpolationCoefficients(uqty);
}

/**
 * Rebuild the interpolation coefficients.
 */
void AdvectionDiffusionTerm::RebuildInterpolationCoefficients(UnknownQuantityHandler* uqty) {
    deltar->SetCoefficient(this->fr, this->drr, uqty, advectionInterpolationMethod_r,  dampingWithIteration*fluxLimiterDampingFactor);
    delta1->SetCoefficient(this->f1, this->d11, uqty, advectionInterpolationMethod_p1, dampingWithIteration*fluxLimiterDampingFactor);
    delta2->SetCoefficient(this->f2, this->d22, uqty, advectionInterpolationMethod_p2, dampingWithIteration*fluxLimiterDampingFactor);
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

