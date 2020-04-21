/**
 * Implementation of a general diffusion term.
 */

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * rg:                Grid on which to define this diffusion term.
 * allocCoefficients: If 'true', allocates memory for the diffusion
 *                    coefficients owned by this term. If 'false',
 *                    it is expected that the memory to use for the
 *                    diffusion coefficients is set using a call to
 *                    'SetCoefficients()' immediately after creation.
 */
DiffusionTerm::DiffusionTerm(Grid *rg, bool allocCoefficients)
    : EquationTerm(rg) {

    if (allocCoefficients)
        this->AllocateCoefficients();
}

/**
 * Destructor.
 */
DiffusionTerm::~DiffusionTerm() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();
}

/**
 * Allocate new memory for the diffusion coefficients.
 */
void DiffusionTerm::AllocateCoefficients() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();
    
    this->drr = new real_t*[nr+1];
    this->d11 = new real_t*[nr];
    this->d12 = new real_t*[nr];
    this->d21 = new real_t*[nr];
    this->d22 = new real_t*[nr];

    len_t
        nElements_fr = n1[nr-1]*n2[nr-1],
        nElements_f1 = 0,
        nElements_f2 = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_fr += n1[i]*n2[i];
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
    }

    this->drr[0] = new real_t[nElements_fr];
    this->d11[0] = new real_t[nElements_f1];
    this->d12[0] = new real_t[nElements_f1];
    this->d22[0] = new real_t[nElements_f2];
    this->d21[0] = new real_t[nElements_f2];

    for (len_t i = 1; i < nr; i++) {
        this->drr[i] = this->drr[i-1] + (n1[i-1]*n2[i-1]);
        this->d11[i] = this->d11[i-1] + ((n1[i-1]+1)*n2[i-1]);
        this->d12[i] = this->d12[i-1] + ((n1[i-1]+1)*n2[i-1]);
        this->d22[i] = this->d22[i-1] + (n1[i-1]*(n2[i-1]+1));
        this->d21[i] = this->d21[i-1] + (n1[i-1]*(n2[i-1]+1));
    }

    // XXX Here we explicitly assume that n1[i] = n1[i+1]
    // at all radii
    //this->drr[nr] = new real_t[n1[nr-1]*n2[nr-1]];
    this->drr[nr] = this->drr[nr-1] + (n1[nr-1]*n2[nr-1]);

    this->coefficientsShared = false;
}

/**
 * Deallocates the memory used by the diffusion coefficients.
 */
void DiffusionTerm::DeallocateCoefficients() {
    if (drr != nullptr) {
        for (len_t i = 0; i < grid->GetNr()+1; i++)
            delete [] drr[i];
        delete [] drr;
    }
    if (d11 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d11[i];
        delete [] d11;
    }
    if (d12 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d12[i];
        delete [] d12;
    }
    if (d21 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d21[i];
        delete [] d21;
    }
    if (d22 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d22[i];
        delete [] d22;
    }
}

/**
 * Assign the memory regions to store the coefficients
 * of this term. This means that we will assume that the
 * memory region is 'shared', and will leave it to someone
 * else to 'free' the memory later on.
 *
 * drr: List of R/R diffusion coefficients.
 * d11: List of p1/p1 diffusion coefficients.
 * d12: List of p1/p2 diffusion coefficients.
 * d21: List of p2/p1 diffusion coefficients.
 * d22: List of p2/p2 diffusion coefficients.
 */
void DiffusionTerm::SetCoefficients(
    real_t **drr,
    real_t **d11, real_t **d12,
    real_t **d21, real_t **d22
) {
    DeallocateCoefficients();

    this->drr = drr;
    this->d11 = d11;
    this->d12 = d12;
    this->d21 = d21;
    this->d22 = d22;

    this->coefficientsShared = true;
}

/**
 * This function is called whenever the computational grid is
 * re-built, in case the grid has been re-sized (in which case
 * we might need to re-allocate memory for the diffusion coefficients)
 */
bool DiffusionTerm::GridRebuilt() {
    this->EquationTerm::GridRebuilt();

    // Do not re-build if our coefficients are owned by someone else
    if (this->coefficientsShared)
        return false;

    this->AllocateCoefficients();

    return true;
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 * NOTE: This routine assumes that the diffusion coefficients
 * are independent of all other unknown quantities (solved
 * for at the same time).
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 */
void DiffusionTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *mat
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(mat, nullptr);
}

/**
 * Build the matrix elements for this operator.
 *
 * mat: Matrix to build elements of.
 * rhs: Right-hand-side of equation (not side).
 */
void DiffusionTerm::SetMatrixElements(Matrix *mat, real_t*) {
    #define f(K,I,J,V) mat->SetElement(offset+j*np1+i, offset + ((K)-ir)*np2*np1 + (J)*np1 + (I), (V))
    #   include "DiffusionTerm.set.cpp"
    #undef f
}

/**
 * Instead of building a linear operator (matrix) to apply to a vector
 * 'x', this routine builds immediately the resulting vector.
 *
 * vec: Vector to set elements of.
 * x:   Input x vector.
 */
void DiffusionTerm::SetVectorElements(real_t *vec, const real_t *x) {
    #define f(K,I,J,V) vec[offset+j*np1+i] += (V)*x[offset+((K)-ir)*np2*np1 + (J)*np1 + (I)]
    #   include "DiffusionTerm.set.cpp"
    #undef f
}

/**
 * Save a list of the DiffusionTerm coefficients to a file,
 * specified by the given SFile object.
 */
void DiffusionTerm::SaveCoefficientsSFile(const std::string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    SaveCoefficientsSFile(sf);
    sf->Close();
}
void DiffusionTerm::SaveCoefficientsSFile(SFile *sf) {
    sfilesize_t dims[3];
    // XXX here we assume that all momentum grids are the same
    const sfilesize_t
        nr = this->grid->GetNr(),
        n1 = this->grid->GetMomentumGrid(0)->GetNp1(),
        n2 = this->grid->GetMomentumGrid(0)->GetNp2();

    dims[0]=nr+1; dims[1]=n2; dims[2]=n1;
    if (this->drr != nullptr)
        sf->WriteMultiArray("Drr", this->drr[0], 3, dims);

    dims[0]=nr; dims[1]=n2+1; dims[2]=n1;
    if (this->d21 != nullptr)
        sf->WriteMultiArray("D21", this->d21[0], 3, dims);
    if (this->d22 != nullptr)
        sf->WriteMultiArray("D22", this->d21[0], 3, dims);

    dims[0]=nr; dims[1]=n2; dims[2]=n1+1;
    if (this->d12 != nullptr)
        sf->WriteMultiArray("D12", this->d12[0], 3, dims);
    if (this->d11 != nullptr)
        sf->WriteMultiArray("D11", this->d11[0], 3, dims);
}

