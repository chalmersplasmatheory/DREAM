/**
 * Implementation of a general advection term.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * rg:          Grid on which this term should be defined.
 * allocCoeffs: If 'true', allocates memory for the coefficients
 *              held by this term. Otherwise, it is expected that
 *              coefficient arrays will be given to this term via
 *              calls to 'SetCoefficients()' and
 *              'SetInterpolationCoefficients()' immediately after
 *              creation.
 */
AdvectionTerm::AdvectionTerm(Grid *rg, bool allocCoeffs)
    : EquationTerm(rg) {
    
    if (allocCoeffs) {
        this->AllocateCoefficients();
        this->AllocateInterpolationCoefficients();
    }
}

/**
 * Destructor.
 */
AdvectionTerm::~AdvectionTerm() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();

    if (!this->interpolationCoefficientsShared)
        DeallocateInterpolationCoefficients();
}

/**
 * Allocate new memory for the advection coefficients
 * based on the current grid sizes.
 */
void AdvectionTerm::AllocateCoefficients() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();

    this->fr = new real_t*[nr+1];
    this->f1 = new real_t*[nr];
    this->f2 = new real_t*[nr];

    len_t
        nElements_fr = n1[nr-1]*n2[nr-1],
        nElements_f1 = 0,
        nElements_f2 = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_fr += n1[i]*n2[i];
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
    }

    this->fr[0] = new real_t[nElements_fr];
    this->f1[0] = new real_t[nElements_f1];
    this->f2[0] = new real_t[nElements_f2];

    for (len_t i = 1; i < nr; i++) {
        this->fr[i] = this->fr[i-1] + (n1[i-1]*n2[i-1]);
        this->f1[i] = this->f1[i-1] + ((n1[i-1]+1)*n2[i-1]);
        this->f2[i] = this->f2[i-1] + (n1[i-1]*(n2[i-1]+1));
    }

    // TODO What about this point???
    //this->fr[nr] = new real_t[???];
    // XXX: Here we assume that the momentum grid is the same
    // at all radial grid points, so that n1_{nr+1/2} = n1_{nr-1/2}
    // (and the same for n2)
    this->fr[nr] = this->fr[nr-1];

    this->ResetCoefficients();

    this->coefficientsShared = false;
}

/**
 * Allocate new memory for the interpolation coefficients.
 */
void AdvectionTerm::AllocateInterpolationCoefficients() {
    if (!this->interpolationCoefficientsShared)
        DeallocateInterpolationCoefficients();

    this->deltar = new real_t*[nr];
    this->delta1 = new real_t*[nr];
    this->delta2 = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        const len_t N = n1[i]*n2[i];

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
}

/**
 * Deallocates the memory used by the advection coefficients.
 */
void AdvectionTerm::DeallocateCoefficients() {
    if (f2 != nullptr) {
        delete [] f2[0];
        delete [] f2;
    }
    if (f1 != nullptr) {
        delete [] f1[0];
        delete [] f1;
    }
    if (fr != nullptr) {
        delete [] fr[0];
        delete [] fr;
    }
}

/**
 * Deallocate the memory used by the interpolation coefficients.
 */
void AdvectionTerm::DeallocateInterpolationCoefficients() {
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
 * Assign the memory regions to store the coefficients
 * of this term. This means that we will assume that the
 * memory region is 'shared', and will leave it to someone
 * else to 'free' the memory later on.
 *
 * fr: List of radial advection coefficients.
 * f1: List of first momentum advection coefficients.
 * f2: List of second momentum advection coefficients.
 */
void AdvectionTerm::SetCoefficients(real_t **fr, real_t **f1, real_t **f2) {
    DeallocateCoefficients();

    this->fr = fr;
    this->f1 = f1;
    this->f2 = f2;

    this->coefficientsShared = true;
}

/**
 * Set the interpolation coefficients explicitly.
 * This equation term will then rely on the owner of these
 * coefficients to de-allocate them later on.
 */
void AdvectionTerm::SetInterpolationCoefficients(
    real_t **dr, real_t **d1, real_t **d2
) {
    DeallocateInterpolationCoefficients();

    this->deltar = dr;
    this->delta1 = d1;
    this->delta2 = d2;

    this->interpolationCoefficientsShared = true;
}

/**
 * This function is called whenever the computational grid is
 * re-built, in case the grid has been re-sized (in which case we might
 * need to re-allocate memory for the advection coefficients)
 */
bool AdvectionTerm::GridRebuilt() {
    bool rebuilt = false;
    this->EquationTerm::GridRebuilt();

    // Do not re-build if our coefficients are owned by someone else
    if (!this->coefficientsShared) {
        this->AllocateCoefficients();
        rebuilt = true;
    }
    if (!this->interpolationCoefficientsShared) {
        this->AllocateInterpolationCoefficients();
        rebuilt = true;
    }
    
    return rebuilt;
}

/**
 * Set all advection coefficients to zero.
 */
void AdvectionTerm::ResetCoefficients() {
    const len_t
        nr = this->grid->GetNr();

    for (len_t ir = 0; ir < nr+1; ir++) {
        // XXX here we assume that all momentum grids are the same
        const len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();

        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1; i++)
                this->fr[ir][j*np1 + i] = 0;
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1+1; i++)
                this->f1[ir][j*(np1+1) + i] = 0;
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2+1; j++)
            for (len_t i = 0; i < np1; i++)
                this->f2[ir][j*np1 + i] = 0;
    }
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 * NOTE: This routine assumes that the advection coefficients
 * are independent of all other unknown quantities (solved
 * for at the same time).
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 */
void AdvectionTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *mat
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(mat, nullptr);
}

/**
 * Build the matrix elements for this operator.
 *
 * mat: Matrix to build elements of.
 * rhs: Right-hand-side of equation (not used).
 */
void AdvectionTerm::SetMatrixElements(Matrix *mat, real_t*) {
    #define f(K,I,J,V) mat->SetElement(offset+j*np1+i, offset + ((K)-ir)*np2*np1 + (J)*np1 + (I), (V))
    #   include "AdvectionTerm.set.cpp"
    #undef f
}

/**
 * Instead of building a linear operator (matrix) to apply to a vector
 * 'x', this routine builds immediately the resulting vector.
 *
 * vec: Vector to set elements of.
 * x:   Input x vector.
 */
void AdvectionTerm::SetVectorElements(real_t *vec, const real_t *x) {
    #define f(K,I,J,V) vec[offset+j*np1+i] += (V)*x[offset+((K)-ir)*np2*np1 + (J)*np1 + (I)]
    #   include "AdvectionTerm.set.cpp"
    #undef f
}

/**
 * Save a list of the AdvectionTerm coefficients to a file,
 * specified by the given SFile object.
 */
void AdvectionTerm::SaveCoefficientsSFile(const std::string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    SaveCoefficientsSFile(sf);
    sf->Close();
}
void AdvectionTerm::SaveCoefficientsSFile(SFile *sf) {
    sfilesize_t dims[3];
    // XXX here we assume that all momentum grids are the same
    const sfilesize_t
        nr = this->grid->GetNr(),
        n1 = this->grid->GetMomentumGrid(0)->GetNp1(),
        n2 = this->grid->GetMomentumGrid(0)->GetNp2();

    if (this->fr != nullptr) {
        dims[0]=nr+1; dims[1]=n2; dims[2]=n1;
        sf->WriteMultiArray("Fr", this->fr[0], 3, dims);
    }
    if (this->f2 != nullptr) {
        dims[0]=nr; dims[1]=n2+1; dims[2]=n1;
        sf->WriteMultiArray("F2", this->f2[0], 3, dims);
    }
    if (this->f1 != nullptr) {
        dims[0]=nr; dims[1]=n2; dims[2]=n1+1;
        sf->WriteMultiArray("F1", this->f1[0], 3, dims);
    }
}

