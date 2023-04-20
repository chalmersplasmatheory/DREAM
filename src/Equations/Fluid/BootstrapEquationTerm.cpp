/**
 * Implementaion of base class from which equation terms related to the bootstrap
 * current should be derived.
 */

#include "DREAM/Equations/Fluid/BootstrapEquationTerm.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
BootstrapEquationTerm::BootstrapEquationTerm(
    FVM::Grid* g, FVM::UnknownQuantityHandler* u, IonHandler *ih,
    BootstrapCurrent *bs, OptionConstants::eqterm_bootstrap_bc bc, real_t sf
) : EquationTerm(g), bc(bc), scaleFactor(sf), bs(bs) {

    nZ = ih->GetNZ();
    nzs = ih->GetNzs();

    id_ncold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ions  = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Ni    = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);   // cc summed densities
    id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);
    if (u->HasUnknown(OptionConstants::UQTY_WI_ENER))
        id_Wi = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    AllocateDeltaX();

    isInitialized = false;
}

/**
 * Destructor.
 */
BootstrapEquationTerm::~BootstrapEquationTerm() {
    DeallocateDeltaX();
}

/**
 * Allocate memory for the central differences X_{k+1} - X_{k-1}.
 */
void BootstrapEquationTerm::AllocateDeltaX() {
    if ( (id_X == id_Ni) || (id_X == id_Wi) )
        deltaX = new real_t[nr * nZ];
    else
        deltaX = new real_t[nr];
}

/**
 * Free memory for the central differences.
 */
void BootstrapEquationTerm::DeallocateDeltaX() {
    delete [] deltaX;
}


void BootstrapEquationTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *u) {
    if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_BACKWARDS)
        rebuildBackwardsBC(u);
    else if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_ZERO)
        rebuildZeroBC(u);
}
/**
 * Calculate the central differences.
 * This one uses a backward finite difference for the grid point by the plasma edge.
 */
void BootstrapEquationTerm::rebuildBackwardsBC(FVM::UnknownQuantityHandler *u) {
    real_t *x = u->GetUnknownData(id_X);
    for (len_t ir = 0; ir < nr; ir++)
        if ( (id_X == id_Ni) || (id_X == id_Wi) )
            for (len_t i = ir; i < nr * nZ; i += nr)
                deltaX[i] = ( (ir == nr - 1) ? x[i] : x[i + 1] ) - ( (ir == 0) ? x[i] : x[i - 1] );
        else
            deltaX[ir] = ( (ir == nr-1) ? x[ir]  : x[ir + 1] ) - ( (ir == 0) ? x[ir] : x[ir - 1] );
}

/**
 * Calculate the central differences.
 * This one assumes plasma quantities are zero outside of the edge.
 */
void BootstrapEquationTerm::rebuildZeroBC(FVM::UnknownQuantityHandler *u) {
    real_t *x = u->GetUnknownData(id_X);
    for (len_t ir = 0; ir < nr; ir++)
        if ( (id_X == id_Ni) || (id_X == id_Wi) )
            for (len_t i = ir; i < nr * nZ; i += nr)
                deltaX[i] = ( (ir == nr - 1) ? 0 : x[i + 1] ) - ( (ir == 0) ? x[i] : x[i - 1] );
        else
            deltaX[ir] = ( (ir == nr-1) ? 0 : x[ir + 1] ) - ( (ir == 0) ? x[ir] : x[ir - 1] );
}


/**
 * Sets the Jacobian matrix for the specified block in the given matrix.
 */
bool BootstrapEquationTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*) {
    if (!HasJacobianContribution(derivId))
        return false;
    for (len_t ir = 0; ir < nr; ir++) {
        if (derivId == id_ions)
            for (len_t izs = 0; izs < nzs; izs++) {
                real_t dC = GetPartialCoefficient(ir, derivId, izs, 0);
                setJacobianElement(derivId, jac, ir, nr * izs + ir, dC * deltaX[ir]);
            }
        else
            if ( (derivId == id_Ni) || (derivId == id_Wi) )
                for (len_t iZ = 0; iZ < nZ; iZ++) {
                    real_t dC = GetPartialCoefficient(ir, derivId, 0, iZ);
                    setJacobianElement(derivId, jac, ir, nr * iZ + ir, dC * deltaX[ir]);
                }
            else {
                real_t dC = GetPartialCoefficient(ir, derivId, 0, 0);
                setJacobianElement(derivId, jac, ir, ir, dC * deltaX[ir]);
            }
    }
    return true;
}
void BootstrapEquationTerm::setJacobianElement(len_t derivId, FVM::Matrix *jac, len_t ir, len_t jr, real_t diagonal) {
    if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_BACKWARDS)
        setJacobianElementBackwardsBC(derivId, jac, ir, jr, diagonal);
    else if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_ZERO)
        setJacobianElementZeroBC(derivId, jac, ir, jr, diagonal);
}
/**
 * Helper function.
 * Sets a single Jacobian matrix element.
 * This one uses a backward finite difference for the grid point by the plasma edge.
 */
void BootstrapEquationTerm::setJacobianElementBackwardsBC(len_t derivId, FVM::Matrix *jac, len_t ir, len_t jr, real_t diagonal) {
    if (derivId == id_X) {
        real_t offDiagonal = GetCoefficient(ir, 0);
        if (ir == 0)
            diagonal -= offDiagonal;
        else
            jac->SetElement(ir, jr - 1, -scaleFactor * offDiagonal);
        if (ir == nr - 1)
            diagonal += offDiagonal;
        else
            jac->SetElement(ir, jr + 1, scaleFactor * offDiagonal);
    }
    jac->SetElement(ir, jr, scaleFactor * diagonal);
}
/**
 * Helper function.
 * Sets a single Jacobian matrix element.
 * This one assumes plasma quantities are zero outside of the edge.
 */
void BootstrapEquationTerm::setJacobianElementZeroBC(len_t derivId, FVM::Matrix *jac, len_t ir, len_t jr, real_t diagonal) {
    if (derivId == id_X) {
        real_t offDiagonal = GetCoefficient(ir, 0);
        if (ir == 0)
            diagonal -= offDiagonal;
        else
            jac->SetElement(ir, jr - 1, -scaleFactor * offDiagonal);
        if (ir != nr-1)
            jac->SetElement(ir, jr + 1, scaleFactor * offDiagonal);
    }
    jac->SetElement(ir, jr, scaleFactor * diagonal);
}



/**
 * Set the linear operator matrix elements corresponding to this term.
 */
void BootstrapEquationTerm::SetMatrixElements(FVM::Matrix *mat, real_t* /* rhs */) {
    for (len_t ir = 0; ir < nr; ir++) {
        if ( (id_X == id_Ni) || (id_X == id_Wi) ) {
            for (len_t iZ = 0; iZ < nZ; iZ++)
                setMatrixElement(mat, nr * iZ + ir, GetCoefficient(ir, iZ));
        } else
            setMatrixElement(mat, ir, GetCoefficient(ir, 0));
    }
}
void BootstrapEquationTerm::setMatrixElement(FVM::Matrix *mat, len_t ir, real_t weight) {
    if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_BACKWARDS)
        setMatrixElementBackwardsBC(mat, ir, weight);
    else if (bc == OptionConstants::EQTERM_BOOTSTRAP_BC_ZERO)
        setMatrixElementZeroBC(mat, ir, weight);
}
/**
 * Helper function.
 * Sets a single linear operator matrix element.
 * This one uses a backward finite difference for the grid point by the plasma edge.
 */
void BootstrapEquationTerm::setMatrixElementBackwardsBC(FVM::Matrix *mat, len_t ir, real_t weight) {
    if (ir == 0)
        mat->SetElement(ir, ir, -scaleFactor * weight);
    else
        mat->SetElement(ir, ir - 1, -scaleFactor * weight);
    if (ir == nr - 1)
        mat->SetElement(ir, ir, scaleFactor * weight);
    else
        mat->SetElement(ir, ir + 1, scaleFactor * weight);
}
/**
 * Helper function.
 * Sets a single linear operator matrix element.
 * This one assumes plasma quantities are zero outside of the edge.
 */
void BootstrapEquationTerm::setMatrixElementZeroBC(FVM::Matrix *mat, len_t ir, real_t weight) {
    if (ir == 0)
        mat->SetElement(ir, ir, -scaleFactor * weight);
    else
        mat->SetElement(ir, ir - 1, -scaleFactor * weight);
    if (ir != nr - 1)
        mat->SetElement(ir, ir + 1, scaleFactor * weight);
}



/**
 * Set function vector for this term. For the nonlinear solver.
 */
void BootstrapEquationTerm::SetVectorElements(real_t *vec, const real_t* /* dt */) {
    if (!isInitialized) {
        bs->Rebuild();
        isInitialized = true;
    }
    for (len_t ir = 0; ir < nr; ir++) {
        if ( (id_X == id_Ni) || (id_X == id_Wi) ) {
            for (len_t iZ = 0; iZ < nZ; iZ++)
                vec[ir] += scaleFactor * GetCoefficient(ir, iZ) * deltaX[nr * iZ + ir];
        } else
            vec[ir] += scaleFactor * GetCoefficient(ir, 0) * deltaX[ir];
    }
}
