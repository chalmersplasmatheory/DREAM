/**
 * Evaluates norms of vectors representing solutions to a
 * system of equations. The vectors only contain data for
 * the non-trivial unknown quantities, and norms should be
 * reported for each non-trivial unknown separately.
 */

#include <vector>
#include "FVM/NormEvaluator.hpp"

using namespace std;
using namespace DREAM::FVM;

/**
 * Constructor.
 */
NormEvaluator::NormEvaluator(UnknownQuantityHandler *uqh, const vector<len_t> &nontrivials)
    : unknowns(uqh), nontrivials(nontrivials) {

    this->nNontrivials = nontrivials.size();
    this->retvec = new real_t[this->nNontrivials];
}

/**
 * Destructor.
 */
NormEvaluator::~NormEvaluator() {
    delete [] this->retvec;
}

/**
 * Evaluate the 2-norm for the given vector.
 */
const real_t *NormEvaluator::Norm2(const real_t *vec) {
    this->Norm2(vec, this->retvec);
    return this->retvec;
}
void NormEvaluator:: Norm2(const real_t *vec, real_t *rvec) {
    len_t offset = 0, i = 0;
    for (auto id : this->nontrivials) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(id);
        const len_t N = uqn->NumberOfElements();

        rvec[i] = 0;
        for (len_t j = 0; j < N; j++)
            rvec[i] += vec[offset+j]*vec[offset+j];

        rvec[i] = sqrt(rvec[i]);

        offset += N;
        i++;
    }
}

/**
 * Evaluates the 2-norm of the difference of the two
 * given vectors:
 *
 *   || v1-v2 ||  =  sqrt[ sum_i[ (v1i - v2i)^2 ] ]
 */
const real_t *NormEvaluator::Norm2Diff(const real_t *vec1, const real_t *vec2) {
    this->Norm2Diff(vec1, vec2, this->retvec);
    return this->retvec;
}
void NormEvaluator::Norm2Diff(const real_t *vec1, const real_t *vec2, real_t *rvec) {
    len_t offset = 0, i = 0;
    for (auto id : this->nontrivials) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(id);
        const len_t N = uqn->NumberOfElements();

        rvec[i] = 0;
        for (len_t j = 0; j < N; j++) {
            real_t v = vec1[offset+j] - vec2[offset+j];
            rvec[i] += v*v;
        }

        rvec[i] = sqrt(rvec[i]);

        offset += N;
        i++;
    }
}

