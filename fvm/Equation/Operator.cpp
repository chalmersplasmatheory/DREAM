/**
 * Implementation of the 'Operator' class, which represents
 * a set of EquationTerms operating on the same quantity.
 */

#include <algorithm>
#include "FVM/Equation/Operator.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 */
Operator::Operator(Grid *grid) : grid(grid) {}

/**
 * Destructor.
 */
Operator::~Operator() {}


/**
 * Add an advection term to this operator.
 */
void Operator::AddTerm(AdvectionTerm *a) {
    if (adterm == nullptr)
        adterm = new AdvectionDiffusionTerm(this->grid);

    adterm->Add(a);
    CheckConsistency();
}

/**
 * Add a diffusion term to this operator.
 */
void Operator::AddTerm(DiffusionTerm *d) {
    if (adterm == nullptr)
        adterm = new AdvectionDiffusionTerm(this->grid);

    adterm->Add(d);
    CheckConsistency();
}

/**
 * Set the predetermined parameter of this operator.
 */
void Operator::AddTerm(PredeterminedParameter *p) {
    if (predetermined != nullptr)
        throw OperatorException("A predetermined parameter has already been applied to this quantity.");

    predetermined = p;
    CheckConsistency();
}

/**
 * Add an evaluable term to this operator.
 */
void Operator::AddTerm(EvaluableEquationTerm *t)  {
    eval_terms.push_back(t);
    CheckConsistency();
}

/**
 * Add a general equation term to this operator.
 */
void Operator::AddTerm(EquationTerm *t)  {
    terms.push_back(t);
    CheckConsistency();
}

/**
 * Add a boundary condition to this operator.
 */
void Operator::AddBoundaryCondition(BC::BoundaryCondition *bc) {
    boundaryConditions.push_back(bc);
}

/**
 * Check the consistency of the terms included in this operator.
 */
void Operator::CheckConsistency() {
    if (predetermined != nullptr) {
        if (adterm != nullptr || terms.size() > 0 || boundaryConditions.size() > 0 || eval_terms.size() > 0)
            throw OperatorException("A predetermined quantity cannot have other equation terms.");
    }
}

/**
 * Evaluate the terms of this operator.
 *
 * vec: Vector to store value in.
 * x:   Current value of the unknown quantity to which this
 *      operator is applied.
 */
void Operator::Evaluate(real_t *vec, const real_t *x) {
    if (IsPredetermined()) {
        const real_t *data = this->predetermined->GetData();
        for (len_t i = 0; i < this->grid->GetNCells(); i++)
            vec[i] += data[i];
    } else {
        this->SetVectorElements(vec, x);
    }
}

/**
 * If this operator contains a single EvaluableTerm, applies
 * the transform of that EvaluableTerm to the given vector.
 *
 * Essentially, it is assumed that this operator represents
 * a single (invertible) operator f(x), operating on the
 * unknown quantity 'x'. Further, it is assumed that this term
 * is part of an equation of the form
 *
 *   f(x) + g(y) = 0
 *
 * where 'y' represents one or more unknown quantities, NOT
 * INCLUDING x. Given the vector 'vec', which is assumed to
 * contain 'g(y)', this method calculates
 *
 *   x = f^-1( -g(y) ),
 *
 * where 'f^-1' denotes the inverse of 'f(x)'.
 */
void Operator::EvaluableTransform(real_t *vec) {
    if (this->eval_terms.size() != 1)
        throw OperatorException(
            "This operator must have exactly one evaluable term for it to be evaluable."
        );

    eval_terms[0]->EvaluableTransform(vec);
}

/**
 * Returns true if this operator contains no equation terms.
 */
bool Operator::IsEmpty() const {
    return (
        this->adterm == nullptr &&
        this->predetermined == nullptr &&
        this->eval_terms.empty() &&
        this->terms.empty() &&
        this->boundaryConditions.empty()
    );
}

/**
 * Returns true if this operator is evaluable, i.e. if it
 * consists of exactly one EvaluableEquationTerm. This means
 * that we can solve for the unknown quantity to which this
 * operator is applied in the manner described in the comment
 * for 'EvaluableTransform()' above.
 *
 * We also return 'true' for 'PredeterminedParameter's, since
 * these have an implied IdentityTerm, which is evaluable, in them.
 */
bool Operator::IsEvaluable() const {
    if (this->adterm != nullptr ||
        boundaryConditions.size() > 0 ||
        terms.size() > 0)
        return false;
    else return
        (this->predetermined != nullptr && this->eval_terms.size() == 0) ||
        (this->predetermined == nullptr && this->eval_terms.size() == 1);
}

/**
 * Returns the number of non-zero elements inserted into
 * a linear operator matrix by this operator object.
 */
len_t Operator::GetNumberOfNonZerosPerRow() const {
    len_t nnz = 0;

    //if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow());
    if (this->adterm != nullptr) nnz = max(nnz, adterm->GetNumberOfNonZerosPerRow());
    if (this->predetermined != nullptr) nnz = max(nnz, predetermined->GetNumberOfNonZerosPerRow());

    for (EquationTerm *term : terms)
        nnz = max(nnz, term->GetNumberOfNonZerosPerRow());
    for (EvaluableEquationTerm *term : eval_terms)
        nnz = max(nnz, term->GetNumberOfNonZerosPerRow());

    // Boundary conditions
    for (BC::BoundaryCondition *bc : boundaryConditions)
        nnz += bc->GetNumberOfNonZerosPerRow();
    
    return nnz;
}

/**
 * Returns the number of non-zero elements inserted into
 * a jacobian matrix by this operator object.
 */
len_t Operator::GetNumberOfNonZerosPerRow_jac() const {
    len_t nnz = 0;

    //if (this->tterm != nullptr) nnz = max(nnz, tterm->GetNumberOfNonZerosPerRow_jac());
    if (this->adterm != nullptr) nnz += adterm->GetNumberOfNonZerosPerRow_jac();
    if (this->predetermined != nullptr) nnz += predetermined->GetNumberOfNonZerosPerRow_jac();

    for (auto it = terms.begin(); it != terms.end(); it++)
        nnz += (*it)->GetNumberOfNonZerosPerRow_jac();
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
        nnz += (*it)->GetNumberOfNonZerosPerRow_jac();

    // Boundary conditions
    for (BC::BoundaryCondition *bc : boundaryConditions)
        nnz += bc->GetNumberOfNonZerosPerRow_jac();
    
    return nnz;
}

/**
 * Make the given equation term identifiable.
 *
 * id:  ID with which to identify the given EquationTerm.
 * trm: Pointer to EquationTerm which should be identifiable.
 */
void Operator::MakeIdentifiable(int_t id, EquationTerm *trm) {
    if (identifiableTerms.find(id) != identifiableTerms.end())
        throw OperatorException("An identifiable term with ID " INT_T_PRINTF_FMT " has already been specified.", id);

    identifiableTerms[id] = trm;
}

/**
 * Build the coefficients of all terms in this operator.
 *
 * t:  Time for which to build the operator.
 * dt: Length of time step to take.
 */
void Operator::RebuildTerms(const real_t t, const real_t dt, UnknownQuantityHandler *uqty) {
    // Predetermined value
    if (predetermined != nullptr) {
        this->predetermined->Rebuild(t, dt, uqty);
        return;
    }

    // Evaluatable equation terms
    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Other equation terms
    for (auto it = terms.begin(); it != terms.end(); it++)
        (*it)->Rebuild(t, dt, uqty);

    // Advection-diffusion term
    if (adterm != nullptr)
        adterm->Rebuild(t, dt, uqty);

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
        (*it)->Rebuild(t, uqty);
}

/**
 * Set the specified block in the given jacobian matrix.
 *
 * uqtyId:     ID of the unknown quantity to which the matrix row belongs.
 * derivId:    ID of the unknown quantity with respect to which the
 *             operator should be differentiated.
 * jac:        Jacobian matrix (block) to set.
 * x:          Value of the unknown quantity.
 * printTerms: Print info about which terms contribute.
 */
bool Operator::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x, bool
#ifndef NDEBUG
    printTerms
#endif
) {
    bool contributes = false;

    for (auto it = eval_terms.begin(); it != eval_terms.end(); it++) {
        bool c = (*it)->SetJacobianBlock( uqtyId, derivId, jac, x);
        contributes |= c;
#ifndef NDEBUG
        if (c && printTerms) printf("Contribution from %s", (*it)->GetName().c_str());
#endif
    }

    for (auto it = terms.begin(); it != terms.end(); it++) {
        bool c = (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);
        contributes |= c;
#ifndef NDEBUG
        if (c && printTerms) printf("Contribution from %s", (*it)->GetName().c_str());
#endif
    }

    // Advection-diffusion term?
    if (adterm != nullptr) {
        contributes |=
#ifndef NDEBUG
            adterm->SetJacobianBlock(uqtyId, derivId, jac, x, printTerms);
#else
            adterm->SetJacobianBlock(uqtyId, derivId, jac, x);
#endif
    }

    // Boundary conditions
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++) {
            bool c = (*it)->AddToJacobianBlock(uqtyId, derivId, jac, x);
            contributes |= c;
#ifndef NDEBUG
        if (c && printTerms) printf("Contribution from %s", (*it)->GetName().c_str());
#endif
    }

    return contributes;
}

/**
 * Routine specifically designed for boundary conditions which need to
 * overwrite elements in the jacobian matrix.
 *
 * uqtyId:  ID of the unknown quantity to which the matrix row belongs.
 * derivId: ID of the unknown quantity with respect to which the
 *          operator should be differentiated.
 * jac:     Jacobian matrix (block) to set.
 * x:       Value of the unknown quantity.
 */
bool Operator::SetJacobianBlockBC(
    const len_t uqtyId, const len_t derivId, Matrix *jac, const real_t *x, bool
#ifndef NDEBUG
    printTerms
#endif
) {
    bool contributes = false;
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++) {
        bool c = (*it)->SetJacobianBlock(uqtyId, derivId, jac, x);
        contributes |= c;
#ifndef NDEBUG
        if (c && printTerms) printf("Contribution from %s", (*it)->GetName().c_str());
#endif
    }

    return contributes;
}

/**
 * Set the linear operator matrix elements in the given
 * matrix in order to represent this operator.
 *
 * mat: Matrix to set elements of.
 * rhs: Vector representing equation right-hand-side.
 */
void Operator::SetMatrixElements(Matrix *mat, real_t *rhs) {
    if (this->IsPredetermined()) {
        this->predetermined->SetMatrixElements(mat, rhs);
    } else {
        for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
            (*it)->SetMatrixElements(mat, rhs);

        for (auto it = terms.begin(); it != terms.end(); it++)
            (*it)->SetMatrixElements(mat, rhs);

        // Advection-diffusion term?
        if (adterm != nullptr)
            adterm->SetMatrixElements(mat, rhs);

        // Boundary conditions
        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
            (*it)->AddToMatrixElements(mat, rhs);

        // TODO Partially assemble matrix
        mat->PartialAssemble();

        // Hard set boundary conditions
        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
            (*it)->SetMatrixElements(mat, rhs);

        mat->PartialAssemble();
    }
}

/**
 * Evaluate this operator and assign its value to the
 * given function vector.
 * 
 * vec: Function vector to assign evaluated operator to.
 * x:   Value of the unknown to evaluate the function for.
 */
void Operator::SetVectorElements(real_t *vec, const real_t *x) {
    if (this->IsPredetermined()) {
        this->predetermined->SetVectorElements(vec, x);
    } else {
        for (auto it = eval_terms.begin(); it != eval_terms.end(); it++)
            (*it)->SetVectorElements(vec, x);

        for (auto it = terms.begin(); it != terms.end(); it++)
            (*it)->SetVectorElements(vec, x);

        // Advection-diffusion term?
        if (adterm != nullptr)
            adterm->SetVectorElements(vec, x);

        // Boundary conditions
        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
            (*it)->AddToVectorElements(vec, x);

        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); it++)
            (*it)->SetVectorElements(vec, x);
    }
}

