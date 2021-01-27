/**
 * Implementation of a base class for diagonal equation terms 
 * where the coefficients are a complex function of unknowns:
 *      T = w(U) * x,
 * where w are the weights and U represents all unknown quantities. 
 * Evaluation of weights must be implemented in derived classes. 
 * It is essentially a DiagonalLinearTerm, but where extra support
 * is provided for setting the Jacobian via SetDiffWeights.
 */

#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;


/**
 * Constructor.
 */
DiagonalComplexTerm::DiagonalComplexTerm(Grid *g, UnknownQuantityHandler *u, Grid *operandGrid)
        : DiagonalTerm(g), operandGrid(operandGrid) {
    this->unknowns = u;
}

DiagonalComplexTerm::~DiagonalComplexTerm(){
    DeallocateDiffWeights();
}

/**
* SetDiffWeights sets the weight jacobian dw/dUderiv, assumed to be local
* (so that w at phase-space point z only depends on U(z) and not, for example,
* on integrals of U). 
*/
void DiagonalComplexTerm::AddWeightsJacobian(
    const len_t /*uqtyId*/, const len_t derivId, Matrix *jac, const real_t* x
){
    /**
    * Check if derivId is one of the id's that contributes 
    * to this advection coefficient 
    */
    len_t nMultiples=0;
    if(!HasJacobianContribution(derivId, &nMultiples))
        return;
    
    ResetDiffWeights();
    SetDiffWeights(derivId, nMultiples);

    len_t NCells = grid->GetNCells();

    if (this->operandGrid != nullptr && this->operandGrid != this->grid) {
        const len_t nr = this->grid->GetNr();
        len_t offset = 0;

        bool operandFluid = (this->operandGrid->GetNCells() == nr);
        if (!operandFluid && this->grid->GetNCells() != nr)
            throw FVMException("DiagonalComplexTerm: No support for terms operating on off-diagonal kinetic-kinetic blocks.");

        Grid *g = (operandFluid ? this->grid : this->operandGrid);
        NCells  = g->GetNCells();

        for (len_t n=0; n<nMultiples; n++)
            for (len_t ir = 0; ir < nr; ir++) {
                MomentumGrid *mg = g->GetMomentumGrid(ir);
                const len_t N = mg->GetNCells();

                for (len_t i = 0; i < N; i++)
                    if (operandFluid)
                        jac->SetElement(offset+i, n*NCells+ir, diffWeights[n*NCells + offset+i] * x[ir]);
                    else
                        jac->SetElement(ir, n*NCells+offset+i, diffWeights[n*NCells + offset+i] * x[offset+i]);

                offset += N;
            }
    } else {
        for(len_t n=0; n<nMultiples; n++)
            for(len_t i=0; i<NCells; i++)
                jac->SetElement(i, n*NCells+i, diffWeights[n*NCells + i] * x[i] ); 
    }
}

/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void DiagonalComplexTerm::SetMatrixElements(Matrix *mat, real_t*) {
    SetElementsInternal([&mat](const len_t i, const len_t j, const real_t v) {
        mat->SetElement(i, j, v);
    });
}

/**
 * Set function vector for this term.
 */
void DiagonalComplexTerm::SetVectorElements(real_t *vec, const real_t *x) {
    SetElementsInternal([&vec,&x](const len_t i, const len_t j, const real_t v) {
        vec[i] += v * x[j];
    });
}

/**
 * Internal routine for setting matrix/vector elements.
 */
void DiagonalComplexTerm::SetElementsInternal(
    std::function<void(const len_t, const len_t, const real_t)> f
) {
    if (this->operandGrid != nullptr && this->operandGrid != this->grid) {
        const len_t nr = this->grid->GetNr();
        len_t offset = 0;

        bool operandFluid = (this->operandGrid->GetNCells() == nr);
        if (!operandFluid && this->grid->GetNCells() != nr)
            throw FVMException("DiagonalComplexTerm: No support for terms operating on off-diagonal kinetic-kinetic blocks.");

        Grid *g = (operandFluid ? this->grid : this->operandGrid);

        for (len_t ir = 0; ir < nr; ir++) {
            MomentumGrid *mg = g->GetMomentumGrid(ir);
            const len_t N = mg->GetNCells();

            for (len_t i = 0; i < N; i++)
                if (operandFluid)
                    f(offset+i, ir, weights[offset+i]);
                else
                    f(ir, offset+i, weights[ir]);

            offset += N;
        }
    } else {
        const len_t N = this->grid->GetNCells();
        for (len_t i = 0; i < N; i++)
            f(i, i, weights[i]);
    }
}


/**
 * Set all diffweights to 0.
 */
void DiagonalComplexTerm::ResetDiffWeights(){
    len_t nMultiples = GetMaxNumberOfMultiplesJacobian();
    len_t NCells = grid->GetNCells();

    if (this->operandGrid != nullptr && NCells < this->operandGrid->GetNCells())
        NCells = this->operandGrid->GetNCells();

    for(len_t i = 0; i<NCells*nMultiples; i++){
        diffWeights[i] = 0;
    }
}

/**
 * Allocate differentiation coefficients.
 */
void DiagonalComplexTerm::AllocateDiffWeights() {
    DeallocateDiffWeights();
    len_t nMultiples = GetMaxNumberOfMultiplesJacobian();
    len_t NCells = grid->GetNCells();

    if (this->operandGrid != nullptr && NCells < this->operandGrid->GetNCells())
        NCells = this->operandGrid->GetNCells();

    diffWeights = new real_t[NCells*nMultiples];
    ResetDiffWeights();
}

/**
 * Deallocate differentiation coefficients.
 */
void DiagonalComplexTerm::DeallocateDiffWeights() {
    if(diffWeights != nullptr)
        delete [] diffWeights;
}

