
/**
 * Implementation of an equation term that takes a kinetic equation term
 * and produces the corresponding fluid term obtained when integrating 
 * it over momentum (weighted by Vp/VpVol).
 */

#include "DREAM/Equations/Fluid/KineticEquationTermIntegratedOverMomentum.hpp"
#include <petscmat.h>

using namespace DREAM;

/**
 * Constructor.
 * 
 *  fluidGrid: grid on which this equation term lives
 *  kineticGrid: grid of the kinetic quantity to which kineticOperator belongs
 *  kineticOperator: a kinetic equation term (appearing in an (f,f) block of the equationsystem), the integral of which becomes this equation term
 *  id_f: unknown id of the kinetic quantity (f_hot or f_re)
 *  u: the UnknownQuantityHandler of the EquationSystem which contains this equation term
 *  scaleFactor: rescales all values of this equation term (+1.0 to become exactly the momentum integral of kineticOperator).
 */
KineticEquationTermIntegratedOverMomentum::KineticEquationTermIntegratedOverMomentum(
    FVM::Grid *fluidGrid, FVM::Grid *kineticGrid, FVM::Operator *kineticOperator, const len_t id_f, FVM::UnknownQuantityHandler *u, PetscScalar scaleFactor
) : FVM::EquationTerm(fluidGrid), kineticGrid(kineticGrid), kineticOperator(kineticOperator), id_f(id_f), unknowns(u), scaleFactor(scaleFactor) {
    allocateKineticStorage();
}


/**
 * Destructor
 */
KineticEquationTermIntegratedOverMomentum::~KineticEquationTermIntegratedOverMomentum(){
    deallocateKineticStorage();
}


/**
 * Allocate and initialize grid quantities
 */
void KineticEquationTermIntegratedOverMomentum::allocateKineticStorage(){
    NCells = kineticGrid->GetNCells();
    PetscInt nnzPerRow = kineticOperator->GetNumberOfNonZerosPerRow_jac();
    PetscInt maxMGN = 0;
    for(len_t ir=0; ir<nr; ir++){
        PetscInt Ni = kineticGrid->GetNp1(ir)*kineticGrid->GetNp2(ir);
        if(Ni>maxMGN)
            maxMGN = Ni;
    }
    kineticMatrix     = new FVM::Matrix(NCells,NCells,nnzPerRow);
    integrationMatrix = new FVM::Matrix(nr, NCells, maxMGN);
    kineticVector = new PetscScalar[NCells];
    fluidVector   = new PetscScalar[nr];

    idxFluid = new PetscInt[nr];
    for(PetscInt ir=0; ir<(PetscInt)nr; ir++)
        idxFluid[ir] = ir;

    // build the integration matrix which carries out Grid::IntegralMomentum upon matrix multiplication 'from the right'
    len_t offset=0;
    for(len_t ir=0; ir<nr; ir++){
        FVM::MomentumGrid *mg = kineticGrid->GetMomentumGrid(ir);
        const len_t   np1  = mg->GetNp1(),  np2 = mg->GetNp2();
        const real_t *dp1  = mg->GetDp1(), *dp2 = mg->GetDp2();
        const real_t *Vp   = kineticGrid->GetVp(ir);
        const real_t VpVol = kineticGrid->GetVpVol(ir);
        for (len_t j = 0; j < np2; j++) 
            for (len_t i = 0; i < np1; i++) {
                len_t idx = j*np1 + i;
                integrationMatrix->SetElement(ir, offset + idx, Vp[idx]*dp1[i]*dp2[j]/VpVol);
            }
        offset += np1*np2;
    }
    integrationMatrix->Assemble();
}


/**
 * Deallocator
 */
void KineticEquationTermIntegratedOverMomentum::deallocateKineticStorage(){
    delete kineticMatrix;
    delete integrationMatrix;
    delete [] kineticVector;
    delete [] fluidVector;
    delete [] idxFluid;
}


/**
 * Called when grid is rebuilt. Reallocates memory.
 */
bool KineticEquationTermIntegratedOverMomentum::GridRebuilt(){
    FVM::EquationTerm::GridRebuilt();

    deallocateKineticStorage();
    allocateKineticStorage();

    return true;
}


/**
 * Sets all elements of the helper vectors to zero
 */
void KineticEquationTermIntegratedOverMomentum::zeroVecs(){
    for(len_t i=0; i<NCells; i++)
        kineticVector[i] = 0;
    for(len_t ir=0; ir<nr; ir++)
        fluidVector[ir] = 0;
}


/**
 * Sets matrix elements of this equation term
 */
void KineticEquationTermIntegratedOverMomentum::SetMatrixElements(FVM::Matrix *mat, real_t *rhs){
    // set our temporary storage to zero
    kineticMatrix->Zero();
    zeroVecs(); 

    kineticOperator->SetMatrixElements(kineticMatrix, kineticVector);
    
    // integrate and set rhs
    kineticGrid->IntegralMomentum(kineticVector, fluidVector);
    for(len_t ir=0; ir<nr; ir++)
        rhs[ir] += scaleFactor*fluidVector[ir];

    kineticMatrix->Assemble();
    Mat C;  
    MatMatMult(integrationMatrix->mat(), kineticMatrix->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C); // performs matrix multiplication C = A*B

    // sum over columns: for each column index j, integrate over momentum (i.e. the rows)
    for(PetscInt j=0; j<(PetscInt)NCells; j++){
        for(len_t ir=0; ir<nr; ir++)
            fluidVector[ir] = 0;
        MatGetValues(C,nr,idxFluid,1,&j,fluidVector);
        for(len_t i=0; i<nr; i++)
            mat->SetElement(i,j,scaleFactor*fluidVector[i]);
    }
}


/**
 * Set jacobian block of this equation term
 */
void KineticEquationTermIntegratedOverMomentum::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t *){
    kineticMatrix->Zero();
    const real_t *f = unknowns->GetUnknownData(id_f);
    kineticOperator->SetJacobianBlock(id_f, derivId, kineticMatrix, f);
    kineticMatrix->Assemble();
    Mat C;  
    MatMatMult(integrationMatrix->mat(), kineticMatrix->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C); // performs matrix multiplication C = A*B

    // sum over columns: for each column index j, integrate over momentum (i.e. the rows)
    PetscInt N = unknowns->GetUnknown(derivId)->NumberOfElements();
    for(PetscInt j=0; j<N; j++){
        for(len_t ir=0; ir<nr; ir++)
            fluidVector[ir] = 0;
        MatGetValues(C,nr,idxFluid,1,&j,fluidVector);
        for(len_t i=0; i<nr; i++)
            jac->SetElement(i,j,scaleFactor*fluidVector[i]);
    }
}


/**
 * Set vector elements. 
 */
void KineticEquationTermIntegratedOverMomentum::SetVectorElements(real_t *vec, const real_t*){
    zeroVecs();
    const real_t *f = unknowns->GetUnknownData(id_f);
    // get vector elements from kinetic equation term
    kineticOperator->SetVectorElements(kineticVector,f);
    // integrate kinetic vector over momentum
    kineticGrid->IntegralMomentum(kineticVector, fluidVector);
    for(len_t ir=0; ir<nr;ir++)
        vec[ir] += scaleFactor*fluidVector[ir];
}


/**
 * Returns NNZ from the input kinetic equation term, replacing NumberOfNonZerosPerRow by the kinetic NCells.
 * (assumes that nnz_jac in kineticOperator is of the form nnz + ...)
 */
len_t KineticEquationTermIntegratedOverMomentum::GetNumberOfNonZerosPerRow_jac() const {
    return kineticOperator->GetNumberOfNonZerosPerRow_jac() + NCells - kineticOperator->GetNumberOfNonZerosPerRow();
}
