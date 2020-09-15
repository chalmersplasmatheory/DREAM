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
AdvectionTerm::AdvectionTerm(Grid *g, bool allocCoeffs)
    : EquationTerm(g) {
    
    if (allocCoeffs) {
        this->AllocateCoefficients();
        this->AllocateInterpolationCoefficients();
    }

}

void AdvectionTerm::AllocateInterpolationCoefficients(){
    if (!this->interpolationCoefficientsShared)
        DeallocateInterpolationCoefficients();

    deltar = new AdvectionInterpolationCoefficient(grid, FLUXGRIDTYPE_RADIAL, AdvectionInterpolationCoefficient::AD_BC_MIRRORED, AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);
    delta1 = new AdvectionInterpolationCoefficient(grid, FLUXGRIDTYPE_P1, AdvectionInterpolationCoefficient::AD_BC_DIRICHLET, AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);
    delta2 = new AdvectionInterpolationCoefficient(grid, FLUXGRIDTYPE_P2, AdvectionInterpolationCoefficient::AD_BC_MIRRORED, AdvectionInterpolationCoefficient::AD_BC_MIRRORED);

    this->interpolationCoefficientsShared = false;
}

/**
 * Destructor.
 */
AdvectionTerm::~AdvectionTerm() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();

    DeallocateDifferentiationCoefficients();

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
    this->f1pSqAtZero = new real_t*[nr];

    len_t
        nElements_fr = n1[nr-1]*n2[nr-1],
        nElements_f1 = 0,
        nElements_f2 = 0,
        nElements_f1pSq = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_fr += n1[i]*n2[i];
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
        nElements_f1pSq += n2[i];
    }

    this->fr[0] = new real_t[nElements_fr];
    this->f1[0] = new real_t[nElements_f1];
    this->f2[0] = new real_t[nElements_f2];
    this->f1pSqAtZero[0] = new real_t[nElements_f1pSq];

    for (len_t i = 1; i < nr; i++) {
        this->fr[i] = this->fr[i-1] + (n1[i-1]*n2[i-1]);
        this->f1[i] = this->f1[i-1] + ((n1[i-1]+1)*n2[i-1]);
        this->f2[i] = this->f2[i-1] + (n1[i-1]*(n2[i-1]+1));
        this->f1pSqAtZero[i] = this->f1pSqAtZero[i-1] + n2[i-1];
    }

    // XXX: Here we assume that the momentum grid is the same
    // at all radial grid points, so that n1_{nr+1/2} = n1_{nr-1/2}
    // (and the same for n2)
    this->fr[nr]  = this->fr[nr-1]  + n1[nr-1]*n2[nr-1];

    this->ResetCoefficients();

    this->coefficientsShared = false;
}

/**
 * Allocate differentiation coefficients.
 */
void AdvectionTerm::AllocateDifferentiationCoefficients() {
    DeallocateDifferentiationCoefficients();
    len_t nMultiples = MaxNMultiple();

    this->dfr = new real_t*[(nr+1)*nMultiples];
    this->df1 = new real_t*[nr*nMultiples];
    this->df2 = new real_t*[nr*nMultiples];

    this->df1pSqAtZero = new real_t*[nr*nMultiples];
    this->JacobianColumn = new real_t[grid->GetNCells()];

    len_t
        nElements_fr = n1[nr-1]*n2[nr-1],
        nElements_f1 = 0,
        nElements_f2 = 0,
        nElements_f1pSq = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_fr += n1[i]*n2[i];
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
        nElements_f1pSq += n2[i];
    }
    
    for (len_t n = 0; n<nMultiples; n++){
        this->dfr[n*(nr+1)] = new real_t[nElements_fr];
        this->df1[n*nr] = new real_t[nElements_f1];
        this->df2[n*nr] = new real_t[nElements_f2];
        this->df1pSqAtZero[n*nr] = new real_t[nElements_f1pSq];
    }
    for (len_t n = 0; n<nMultiples; n++){
        for (len_t ir = 1; ir < nr; ir++) {
            this->dfr[ir+n*(nr+1)] = this->dfr[ir-1+n*(nr+1)] + (n1[ir-1]*n2[ir-1]);
            this->df1[ir+n*nr] = this->df1[ir-1+n*nr] + ((n1[ir-1]+1)*n2[ir-1]);
            this->df2[ir+n*nr] = this->df2[ir-1+n*nr] + (n1[ir-1]*(n2[ir-1]+1));
            this->df1pSqAtZero[ir+n*nr] = this->df1pSqAtZero[ir-1+n*nr] + n2[ir-1];
        }

        // XXX: Here we assume that the momentum grid is the same
        // at all radial grid points, so that n1_{nr+1/2} = n1_{nr-1/2}
        // (and the same for n2)
        this->dfr[nr+n*(nr+1)] = this->dfr[nr-1+n*(nr+1)] + n1[nr-1]*n2[nr-1];
    }

    this->ResetDifferentiationCoefficients();
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
    if (f1pSqAtZero != nullptr) {
        delete [] f1pSqAtZero[0];
        delete [] f1pSqAtZero;
    }
}

/**
 * Dellocate the memory used by the differentiation coefficients.
 */
void AdvectionTerm::DeallocateDifferentiationCoefficients() {
    // Differentiation coefficients
    if (df2 != nullptr) {
        delete [] df2[0];
        delete [] df2;
    }
    if (df1 != nullptr) {
        delete [] df1[0];
        delete [] df1;
    }
    if (dfr != nullptr) {
        delete [] dfr[0];
        delete [] dfr;
    }
    if (JacobianColumn != nullptr){
        delete [] JacobianColumn;
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
void AdvectionTerm::SetCoefficients(
    real_t **fr, real_t **f1, real_t **f2, real_t **f1pSqAtZero
) {
    DeallocateCoefficients();

    this->fr = fr;
    this->f1 = f1;
    this->f2 = f2;
    this->f1pSqAtZero = f1pSqAtZero;

    this->coefficientsShared = true;
}

/**
 * Set the interpolation coefficients explicitly.
 * This equation term will then rely on the owner of these
 * coefficients to de-allocate them later on.
 */
void AdvectionTerm::SetInterpolationCoefficients(
    AdvectionInterpolationCoefficient *dr, AdvectionInterpolationCoefficient *d1, AdvectionInterpolationCoefficient *d2
) {
    DeallocateInterpolationCoefficients();

    this->deltar = dr;
    this->delta1 = d1;
    this->delta2 = d2;

    this->interpolationCoefficientsShared = true;
}


void AdvectionTerm::DeallocateInterpolationCoefficients(){
    
    if(deltar!=nullptr){
        delete deltar;
        deltar = nullptr;
    }
    if(delta1!=nullptr){
        delete delta1;
        delta1 = nullptr;
    }
    if(delta2!=nullptr){
        delete delta2;
        delta2 = nullptr;
    }
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
        deltar->GridRebuilt();
        delta1->GridRebuilt();
        delta2->GridRebuilt();
        deltar->SetCoefficient(this->fr);
        delta1->SetCoefficient(this->f1);
        delta2->SetCoefficient(this->f2);

        rebuilt = true;
    }

    // TODO: find condition for when to allocate these
    AllocateDifferentiationCoefficients();
    
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
                this->fr[ir][j*np1 + i]  = 0;
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2; j++){
            for (len_t i = 0; i < np1+1; i++)
                this->f1[ir][j*(np1+1) + i]  = 0;
            this->f1pSqAtZero[ir][j] = 0;
        }
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2+1; j++)
            for (len_t i = 0; i < np1; i++)
                this->f2[ir][j*np1 + i]  = 0;
    }
}

/**
 * Set all differentiation coefficients to zero.
 */
void AdvectionTerm::ResetDifferentiationCoefficients() {
    len_t nMultiples = MaxNMultiple();

    const len_t
        nr = this->grid->GetNr();

    for(len_t n=0; n<nMultiples; n++){
        for (len_t ir = 0; ir < nr+1; ir++) {
            // XXX here we assume that all momentum grids are the same
            const len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();
            const len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();

            for (len_t j = 0; j < np2; j++)
                for (len_t i = 0; i < np1; i++)
                    this->dfr[ir+n*(nr+1)][j*np1 + i] = 0;
        }

        for (len_t ir = 0; ir < nr; ir++) {
            const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
            const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

            for (len_t j = 0; j < np2; j++){
                for (len_t i = 0; i < np1+1; i++)
                    this->df1[ir+n*nr][j*(np1+1) + i] = 0;
                this->df1pSqAtZero[ir+n*nr][j] = 0;
            }

            for (len_t j = 0; j < np2+1; j++)
                for (len_t i = 0; i < np1; i++)
                    this->df2[ir+n*nr][j*np1 + i] = 0;
        }
    }
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * jac:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 */
void AdvectionTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x
) {
    interp_mode = AdvectionInterpolationCoefficient::AD_INTERP_MODE_JACOBIAN;
    if ( (uqtyId == derivId) && !this->coefficientsShared)
        this->SetMatrixElements(jac, nullptr);
    
   /**
    * Check if derivId is one of the id's that contribute
    * to this advection coefficient 
    */
    bool hasDerivIdContribution = false;
    len_t nMultiples;
    for(len_t i_deriv = 0; i_deriv < derivIds.size(); i_deriv++){
        if (derivId == derivIds[i_deriv]){
            nMultiples = derivNMultiples[i_deriv];
            hasDerivIdContribution = true;
        }
    }
    if(!hasDerivIdContribution)
        return;
    

    // TODO: allocate differentiation coefficients in a more logical location
    if(df1 == nullptr)
        AllocateDifferentiationCoefficients();

    // Set partial advection coefficients for this advection term 
    SetPartialAdvectionTerm(derivId, nMultiples);

    len_t offset;
    for(len_t n=0; n<nMultiples; n++){
        ResetJacobianColumn();
        SetVectorElements(JacobianColumn, x, dfr+n*(nr+1), df1+n*nr, df2+n*nr, df1pSqAtZero+n*nr);
        offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            for (len_t j = 0; j < n2[ir]; j++) 
                for (len_t i = 0; i < n1[ir]; i++) 
                    jac->SetElement(offset + n1[ir]*j + i, n*nr+ir, JacobianColumn[offset + n1[ir]*j + i]); 

            offset += n1[ir]*n2[ir];
        }
    }
}

void AdvectionTerm::ResetJacobianColumn(){
    len_t offset = 0; 
    for(len_t ir=0; ir<nr; ir++){
        for (len_t j = 0; j < n2[ir]; j++) 
            for (len_t i = 0; i < n1[ir]; i++) 
                JacobianColumn[offset + n1[ir]*j + i] = 0;

        offset += n1[ir]*n2[ir];
    }

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
    this->SetVectorElements(vec, x, this->fr, this->f1, this->f2, this->f1pSqAtZero);
}
void AdvectionTerm::SetVectorElements(
    real_t *vec, const real_t *x,
    const real_t *const* fr, const real_t *const* f1, const real_t *const* f2, const real_t *const* f1pSqAtZero
) {
    interp_mode = AdvectionInterpolationCoefficient::AD_INTERP_MODE_FULL;
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

