/**
 * An external boundary condition which represents a physically
 * motivated boundary condition for the poloidal flux psi_p at
 * the plasma boundary r=a, which is related to the total plasma 
 * current I_p.
 */

#include "DREAM/Equations/Fluid/AmperesLawBoundaryAtRMax.hpp"
#include "DREAM/Constants.hpp"


using namespace DREAM::FVM::BC;




/**
 * Constructor.
 */
AmperesLawBoundaryAtRMax::AmperesLawBoundaryAtRMax(Grid *g, Grid *targetGrid, 
    const Operator *eqn, real_t scaleFactor)
    : BoundaryCondition(g), equation(eqn), targetGrid(targetGrid), scaleFactor(scaleFactor) {
    
    SetName("AmperesLawBoundaryAtRMax");
}

/**
 * Destructor.
 */
AmperesLawBoundaryAtRMax::~AmperesLawBoundaryAtRMax(){}

/**
 * Rebuilds the single non-zero matrix element that contributes, which 
 * represents the external flux at the upper r boundary in the 
 * AmperesLawDiffusionTerm induced by the total plasma current I_p.
 */
bool AmperesLawBoundaryAtRMax::Rebuild(const real_t, UnknownQuantityHandler*){
    const len_t nr = this->grid->GetNr();    
    RadialGrid *rGrid = grid->GetRadialGrid();
    
    real_t
        dr    = rGrid->GetDr()[nr-1], 
        dr_f  = rGrid->GetR_f(nr) - rGrid->GetR(nr-1), 
        Vp    = this->grid->GetVp(nr-1)[0],
        Vp_fr = this->grid->GetVp_fr(nr)[0],
        Drr = equation->GetDiffusionCoeffRR(nr)[0];

    real_t diffusionTermCoeff = -Drr*Vp_fr/(Vp*dr);

    /**
     *  dpsi/dr(r=a) = [psi(a)-psi(rmax)]/(a-rmax) 
     */
    real_t dPsiDrCoeff = scaleFactor/dr_f;

    this->coefficient = diffusionTermCoeff * dPsiDrCoeff;

    return true;
}

/**
 * Add flux to jacobian block.
 */
bool AmperesLawBoundaryAtRMax::AddToJacobianBlock(
    const len_t derivId, const len_t qtyId, Matrix * jac, const real_t* /*x*/
) {
    if (derivId == qtyId)
        this->AddToMatrixElements(jac, nullptr);
    return (derivId == qtyId);
}

/**
 * Add flux to linearized operator matrix.
 *
 * mat: Matrix to add boundary conditions to.
 * rhs: Right-hand-side vector (not used).
 */
void AmperesLawBoundaryAtRMax::AddToMatrixElements(
    Matrix *mat, real_t*
) {
    const len_t nr = grid->GetNr(); 
    const len_t N_target = targetGrid->GetNCells();   
    mat->SetElement(nr-1, N_target-1, coefficient);
}
/**
 * Add flux to function vector.
 *
 * vec: Function vector to add boundary conditions to.
 * f:   Current value of distribution function.
 */
void AmperesLawBoundaryAtRMax::AddToVectorElements(
    real_t *vec, const real_t *f
) {
    const len_t nr = this->grid->GetNr();
    const len_t N_target = targetGrid->GetNCells();   

    vec[nr-1] += coefficient * f[N_target-1];
}
