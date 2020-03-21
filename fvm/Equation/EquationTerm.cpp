/**
 * Implementation of the 'EquationTerm' base class.
 */

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"

using namespace TQS::FVM;


/**
 * Constructor.
 */
EquationTerm::EquationTerm(RadialGrid *rg, bool allocInterpolationCoeffs)
    : grid(rg) {

    if (allocInterpolationCoeffs)
        this->AllocateInterpolationCoefficients();
}

/**
 * Destructor.
 */
EquationTerm::~EquationTerm() {
    if (!this->interpolationCoeffsShared)
        DeallocateInterpolationCoefficients();
}

/**
 * Allocate memory for the interpolation coefficients.
 */
void EquationTerm::AllocateInterpolationCoefficients() {
    if (!this->interpolationCoeffsShared)
        DeallocateInterpolationCoefficients();

    len_t nr = this->grid->GetNr();

    this->deltar = new real_t*[nr];
    this->delta1 = new real_t*[nr];
    this->delta2 = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        len_t N = this->grid->GetMomentumGrid(i)->GetNCells();

        this->deltar[i] = new real_t[N];
        this->delta1[i] = new real_t[N];
        this->delta2[i] = new real_t[N];

        // Initialize to delta = 1/2
        for (len_t j = 0; j < N; j++) {
            this->deltar[i][j] = 0.5;
            this->delta1[i][j] = 0.5;
            this->delta2[i][j] = 0.5;
        }
    }

    this->interpolationCoeffsShared = false;
}

/**
 * De-allocate memory for the interpolation coefficients.
 */
void EquationTerm::DeallocateInterpolationCoefficients() {
    if (delta2 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] delta2[i];

        delete [] delta2;
    }

    if (delta1 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] delta1[i];

        delete [] delta1;
    }

    if (deltar != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] deltar[i];

        delete [] deltar;
    }
}

/**
 * Function called when any of the grids have been re-built.
 */
bool EquationTerm::GridRebuilt() {
    // Do not re-build if our coefficients are owned by someone else
    if (this->interpolationCoeffsShared)
        return false;

    this->AllocateInterpolationCoefficients();

    return true;
}

/**
 * Set the interpolation coefficients explicitly.
 * This equation term will then rely on the owner of these
 * coefficients to de-allocate them later on.
 */
void EquationTerm::SetInterpolationCoefficients(
    real_t **dr, real_t **d1, real_t **d2
) {
    DeallocateInterpolationCoefficients();

    this->deltar = dr;
    this->delta1 = d1;
    this->delta2 = d2;

    this->interpolationCoeffsShared = true;
}

