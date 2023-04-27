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
    BootstrapCurrent *bs, real_t sf
) : EquationTerm(g), scaleFactor(sf), bs(bs) {

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

    real_t *x = u->GetUnknownData(id_X);

    // set bc
    if ( (id_X == id_Ni) || (id_X == id_Wi) )
        for (len_t i = 0; i < nr * nZ; i += nr) {
            deltaX[i] = x[i+1] - x[i];
            deltaX[i+nr-1] = x[i+nr-1]-x[i+nr-2];
        }
    else {
        deltaX[0] = x[1] - x[0];
        deltaX[nr-1] = x[nr-1] - x[nr-2];
    }

    for (len_t ir = 1; ir < nr - 1; ir++)
        if ( (id_X == id_Ni) || (id_X == id_Wi) )
            for (len_t i = ir; i < nr * nZ; i++)
                deltaX[i] = x[i+1] - x[i-1];
        else
            deltaX[ir] = x[ir+1] - x[ir-1];
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
/**
 * Helper function.
 * Sets a single Jacobian matrix element, with a backward finite difference for the grid point by the plasma edge.
 */
void BootstrapEquationTerm::setJacobianElement(len_t derivId, FVM::Matrix *jac, len_t ir, len_t jr, real_t diagonal) {
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
/**
 * Helper function.
 * Sets a single linear operator matrix element, with a backward finite difference for the grid point by the plasma edge.
 */
void BootstrapEquationTerm::setMatrixElement(FVM::Matrix *mat, len_t ir, real_t weight) {
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
