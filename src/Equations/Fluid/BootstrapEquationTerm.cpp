/**
 * Implementaion of base class from which equation terms related to the bootstrap
 * current should be derived.
 */

using namespace DREAM;

/**
 * Constructor.
 */
BootstrapEquationTerm::BootstrapEquationTerm(
    FVM::Grid* g, FVM::UnknownQuantityHandler* u, IonHandler *ih, real_t sf
) : EquationTerm(g), scaleFactor(sf) {

    nZ = ih->GetNZ();
    nzs = ih->GetNZs();

    id_ncold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ions  = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Ni    = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);   // cc summed densities
    id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);
    if (u->HasUnknown(OptionConstants::UQQTY_WI_ENER))
        id_Wi = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);

    AllocateDeltaX();
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
void AllocateDeltaX() {
    if ( (id_X == id_Ni) || (id_X == id_Wi) )
        deltaX = new real_t[nr * nZ];
    else
        deltaX = new real_t[nr];
}

/**
 * Free memory for the central differences.
 */
void DeallocateDeltaX() {
    delete [] deltaX;
}

/**
 * Calculate the central differences, with boundary conditions:
 *      X_{-1}   = X_{0}
 *      X_{nr-1} = 0
 */
void BootstrapEquationTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* u) {
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
bool BootstrapEquationTerm::SetJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*
) {
    if (!HasJacobianContribution(derivId))
        return false;
    for (len_t ir = 0; ir < nr; ir++) {
        if (derivId == id_ions)
            for (len_t izs = 0; izs < nzs; izs++) {
                real_t dC = GetPartialCoefficient(ir, derivId, izs, nullptr);
                SetJacobianElement(derivId, jac, ir, nr * izs + ir, dC * deltaX[ir]);
            }
        else
            if ( (derivId == id_Ni) || (derivId == id_Wi) )
                for (len_t iZ = 0; iZ < nZ; iZ++) {
                    real_t dC = GetPartialCoefficient(ir, derivId, nullptr, iZ);
                    SetJacobianElement(derivId, jac, ir, nr * iZ + ir, dC * deltaX[ir]);
                }
            else {
                real_t dC = GetPartialCoefficient(ir, derivId, nullptr, nullptr);
                SetJacobianElement(derivId, jac, ir, ir, dC * deltaX[ir]);
            }
    }
    return true;
}

/**
 * Helper function. Sets a single Jacobian matrix element.
 */
void SetJacobianElement(len_t derivId, FVM::Matrix *jac, len_t ir, len_t jr, real_t diagonal) {
    if (derivId == id_X) {
        real_t offDiagonal = GetCoefficient(ir, nullptr);
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
                SetMatrixElement(mat, nr * iZ + ir, GetCoefficient(ir, iZ));
        } else
            SetMatrixElement(mat, ir, GetCoefficient(ir, nullptr))
}

/**
 * Helper function. Sets a single linear operator matrix element.
 */
void BootstrapEquationTerm::SetMatrixElement(FVM::Matrix, len_t ir, real_t weight) {
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
void BootstrapEquationTerm::SetVectorElements(real_t *vec, const real_t*) {
    for (len_t ir = 0; ir < nr; ir++) {
        if ( (id_X == id_Ni) || (id_X == id_Wi) ) {
            for (len_t iZ = 0; iZ < nZ; iZ++)
                vec[ir] += scaleFactor * GetCoefficient(ir, iZ) * deltaX[nr * iZ + ir];
        } else
            vec[ir] += scaleFactor * GetCoefficient(ir, nullptr) * deltaX[ir];
}
