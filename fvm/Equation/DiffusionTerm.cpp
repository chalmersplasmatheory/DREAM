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
    
    DeallocateDifferentiationCoefficients();
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
    this->drr[nr]  = this->drr[nr-1]  + (n1[nr-1]*n2[nr-1]);

    this->ResetCoefficients();

    // Set interpolation coefficients for unknowns to radial flux grid
    this->deltaRadialFlux = new real_t[nr+1];
    for(len_t ir=0; ir<nr+1; ir++){
        if(ir==0)
            deltaRadialFlux[0] = 1;
        else if (ir<nr)
            deltaRadialFlux[ir] = 0.5*grid->GetRadialGrid()->GetDr(ir-1)/grid->GetRadialGrid()->GetDr_f(ir-1); // linear interpolation coefficient
        else 
            deltaRadialFlux[nr] = 0;
    }

    this->coefficientsShared = false;
}

void DiffusionTerm::AllocateDifferentiationCoefficients() {
    DeallocateDifferentiationCoefficients();
    len_t nMultiples = MaxNMultiple();

    this->ddrr = new real_t*[(nr+1)*nMultiples];
    this->dd11 = new real_t*[nr*nMultiples];
    this->dd12 = new real_t*[nr*nMultiples];
    this->dd21 = new real_t*[nr*nMultiples];
    this->dd22 = new real_t*[nr*nMultiples];

    this->JacobianColumn = new real_t[grid->GetNCells()];


    len_t
        nElements_fr = n1[nr-1]*n2[nr-1],
        nElements_f1 = 0,
        nElements_f2 = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_fr += n1[i]*n2[i];
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
    }


    for (len_t n = 0; n<nMultiples; n++){
        this->ddrr[n*(nr+1)] = new real_t[nElements_fr];
        this->dd11[n*nr] = new real_t[nElements_f1];
        this->dd12[n*nr] = new real_t[nElements_f1];
        this->dd22[n*nr] = new real_t[nElements_f2];
        this->dd21[n*nr] = new real_t[nElements_f2];
    }


    for (len_t n = 0; n<nMultiples; n++){
        for (len_t i = 1; i < nr; i++) {
            this->ddrr[i+n*(nr+1)] = this->ddrr[i-1+n*(nr+1)] + (n1[i-1]*n2[i-1]);
            this->dd11[i+n*nr] = this->dd11[i-1+n*nr] + ((n1[i-1]+1)*n2[i-1]);
            this->dd12[i+n*nr] = this->dd12[i-1+n*nr] + ((n1[i-1]+1)*n2[i-1]);
            this->dd22[i+n*nr] = this->dd22[i-1+n*nr] + (n1[i-1]*(n2[i-1]+1));
            this->dd21[i+n*nr] = this->dd21[i-1+n*nr] + (n1[i-1]*(n2[i-1]+1));
        }

        // XXX Here we explicitly assume that n1[i] = n1[i+1]
        // at all radii
        this->ddrr[nr+n*(nr+1)] = this->ddrr[nr-1+n*(nr+1)] + (n1[nr-1]*n2[nr-1]);

    }
    
    this->ResetDifferentiationCoefficients();
}

/**
 * Deallocates the memory used by the diffusion coefficients.
 */
void DiffusionTerm::DeallocateCoefficients() {
    if (drr != nullptr) {
        delete [] drr[0];
        delete [] drr;
    }
    if (d11 != nullptr) {
        delete [] d11[0];
        delete [] d11;
    }
    if (d12 != nullptr) {
        delete [] d12[0];
        delete [] d12;
    }
    if (d21 != nullptr) {
        delete [] d21[0];
        delete [] d21;
    }
    if (d22 != nullptr) {
        delete [] d22[0];
        delete [] d22;
    }

    if(deltaRadialFlux != nullptr)
        delete [] deltaRadialFlux;
}

/**
 * Deallocates the memory used by the diffusion coefficients.
 */
void DiffusionTerm::DeallocateDifferentiationCoefficients() {
    if (ddrr != nullptr) {
        delete [] ddrr[0];
        delete [] ddrr;
    }
    if (dd11 != nullptr) {
        delete [] dd11[0];
        delete [] dd11;
    }
    if (dd12 != nullptr) {
        delete [] dd12[0];
        delete [] dd12;
    }
    if (dd21 != nullptr) {
        delete [] dd21[0];
        delete [] dd21;
    }
    if (dd22 != nullptr) {
        delete [] dd22[0];
        delete [] dd22;
    }

    if(JacobianColumn != nullptr)
        delete [] JacobianColumn;
}

/**
 * Assign the memory regions to store the coefficients
 * of this term. This means that we will assume that the
 * memory region is 'shared', and will leave it to someone
 * else to 'free' the memory later on.
 *
 * drr:  List of R/R diffusion coefficients.
 * d11:  List of p1/p1 diffusion coefficients.
 * d12:  List of p1/p2 diffusion coefficients.
 * d21:  List of p2/p1 diffusion coefficients.
 * d22:  List of p2/p2 diffusion coefficients.
 *
 * ddrr: List of R/R differentiation coefficients.
 * dd11: List of p1/p1 differentiation coefficients.
 * dd12: List of p1/p2 differentiation coefficients.
 * dd21: List of p2/p1 differentiation coefficients.
 * dd22: List of p2/p2 differentiation coefficients.
 */
void DiffusionTerm::SetCoefficients(
    real_t **drr,
    real_t **d11, real_t **d12,
    real_t **d21, real_t **d22,
    real_t *delta
) {
    DeallocateCoefficients();

    this->drr = drr;
    this->d11 = d11;
    this->d12 = d12;
    this->d21 = d21;
    this->d22 = d22;
    this->deltaRadialFlux = delta;

    this->coefficientsShared = true;
}

/**
 * This function is called whenever the computational grid is
 * re-built, in case the grid has been re-sized (in which case
 * we might need to re-allocate memory for the diffusion coefficients)
 */
bool DiffusionTerm::GridRebuilt() {
    this->EquationTerm::GridRebuilt();
    
    // TODO: find condition for when to allocate these
    this->AllocateDifferentiationCoefficients();
    
    // Do not re-build if our coefficients are owned by someone else
    if (this->coefficientsShared)
        return false;

    this->AllocateCoefficients();
    return true;
}

/**
 * Set all advection coefficients to zero.
 */
void DiffusionTerm::ResetCoefficients() {
    const len_t
        nr = this->grid->GetNr();

    for (len_t ir = 0; ir < nr+1; ir++) {
        // XXX here we assume that all momentum grids are the same
        const len_t np2 = this->grid->GetMomentumGrid(0)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(0)->GetNp1();

        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1; i++)
                this->drr[ir][j*np1 + i]  = 0;
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1+1; i++) {
                this->d11[ir][j*(np1+1) + i]  = 0;
                this->d12[ir][j*(np1+1) + i]  = 0;
            }
    }

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

        for (len_t j = 0; j < np2+1; j++)
            for (len_t i = 0; i < np1; i++) {
                this->d22[ir][j*np1 + i]  = 0;
                this->d21[ir][j*np1 + i]  = 0;
            }
    }
}

/**
 * Set all differentiation coefficients to zero.
 */
void DiffusionTerm::ResetDifferentiationCoefficients() {
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
                    this->ddrr[ir+n*(nr+1)][j*np1 + i] = 0;
        }
        
        for (len_t ir = 0; ir < nr; ir++) {
            const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
            const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();

            for (len_t j = 0; j < np2; j++)
                for (len_t i = 0; i < np1+1; i++) {
                    this->dd11[ir+n*nr][j*(np1+1) + i] = 0;
                    this->dd12[ir+n*nr][j*(np1+1) + i] = 0;
                }
        
            for (len_t j = 0; j < np2+1; j++)
                for (len_t i = 0; i < np1; i++) {
                    this->dd22[ir+n*nr][j*np1 + i] = 0;
                    this->dd21[ir+n*nr][j*np1 + i] = 0;
                }
        }
    }
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
 * jac:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 */
void DiffusionTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t* x
) {
    if ( (uqtyId == derivId) && !this->coefficientsShared)
        this->SetMatrixElements(jac, nullptr);

    
    /**
    * Check if derivId is one of the id's that contributes 
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
    if(dd11 == nullptr)
        AllocateDifferentiationCoefficients();

    // Set partial advection coefficients for this advection term 
    SetPartialDiffusionTerm(derivId, nMultiples);

    for(len_t n=0; n<nMultiples; n++){
        SetPartialJacobianContribution(0, JACOBIAN_SET_CENTER, n, jac, x);
        SetPartialJacobianContribution(-1,JACOBIAN_SET_LOWER, n, jac, x);
        SetPartialJacobianContribution(+1,JACOBIAN_SET_UPPER, n, jac, x);
    }
}

/**
 * Sets one of the diagonals in the jacobian block.
 */
void DiffusionTerm::SetPartialJacobianContribution(int_t diagonalOffset, jacobian_interp_mode set_mode, len_t n, Matrix *jac, const real_t *x){
        ResetJacobianColumn();
        SetVectorElements(JacobianColumn, x, ddrr+n*(nr+1),
                            dd11+n*nr, dd12+n*nr,
                            dd21+n*nr, dd22+n*nr, set_mode);
        len_t offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            if((ir==0&&diagonalOffset==-1) || ir+diagonalOffset>=nr)
                continue;
            for (len_t j = 0; j < n2[ir]; j++)
                for (len_t i = 0; i < n1[ir]; i++)
                    jac->SetElement(offset + n1[ir]*j + i, n*nr+ir+diagonalOffset, JacobianColumn[offset + n1[ir]*j + i]);
            offset += n1[ir]*n2[ir];
        }
}

/**
 * Sets the jacobian helper vector to zero
 */
void DiffusionTerm::ResetJacobianColumn(){
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
 * rhs: Right-hand-side of equation (not side).
 */
void DiffusionTerm::SetMatrixElements(Matrix *mat, real_t*) {
    jacobian_interp_mode set = NO_JACOBIAN;
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
    this->SetVectorElements(
        vec, x, this->drr,
        this->d11, this->d12,
        this->d21, this->d22
    );
}
void DiffusionTerm::SetVectorElements(
    real_t *vec, const real_t *x,
    const real_t *const* drr,
    const real_t *const* d11, const real_t *const* d12,
    const real_t *const* d21 ,const real_t *const* d22,
    jacobian_interp_mode set
) {
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

