/**
 * Implementation of the 'Matrix' class.
 * This object represents a view into a part of the full
 * matrix representing the equation system, found in the
 * 'EquationSystem' class.
 */

#include <iostream>
#include <petscmat.h>
#include "FVM/Matrix.hpp"


using namespace std;
using namespace DREAM::FVM;

/**
 * Constructors.
 */
Matrix::Matrix() {}
Matrix::Matrix(const PetscInt m, const PetscInt n, Mat mat) {
    this->m = m;
    this->n = n;
    this->petsc_mat = mat;
}

/**
 * Creates a new Matrix of size m-by-n, with nnz non-zero
 * elements pre-allocated. If 'distributed' is true, the
 * matrix is stored in distributed memory (i.e. across
 * nodes).
 *
 * m:           Number of matrix rows.
 * n:           Number of matrix columns.
 * nnz:         Number of non-zero elements per row expected to be
 *              put in the matrix (to be pre-allocated).
 *              If 'nnzl' is also given, this parameter
 *              only indicates the number of elements in 'nnzl'.
 * nnzl:        List of number of non-zero elements in each
 *              row (can be different for each row). If
 *              'nullptr', 'nnz' is interpreted as the total
 *              number of non-zero elements.
 */
Matrix::Matrix(
    const PetscInt m, const PetscInt n, const PetscInt nnz,
    const PetscInt *nnzl
) { this->Construct(m, n, nnz, nnzl); }

void Matrix::Construct(
    const PetscInt m, const PetscInt n, const PetscInt nnz,
    const PetscInt *nnzl
) {
    PetscErrorCode ierr;

    if (this->allocated)
        this->Destroy();

    this->m = m;
    this->n = n;

    MatCreate(PETSC_COMM_WORLD, &(this->petsc_mat));
    //MatSetType(this->petsc_mat, MATSEQAIJ);
    MatSetType(this->petsc_mat, MATAIJ);
    MatSetSizes(this->petsc_mat, PETSC_DECIDE, PETSC_DECIDE, m, n);

    if ((ierr=MatSeqAIJSetPreallocation(this->petsc_mat, nnz, nnzl)))
        throw MatrixException("Failed to allocate memory for PETSc matrix. Error code: %d", ierr);

    // Ensure that the non-zero structure of the matrix is
    // kept when 'ZeroRows()' is called.
    MatSetOption(this->petsc_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    // Don't complain about new allocations
    MatSetOption(this->petsc_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    this->allocated = true;

    // Set diagonal elements to zero (because PETSc requires all
    // diagonals to be explicitly set, even if not used)
    for (PetscInt i = 0; i < m; i++)
        MatSetValue(this->petsc_mat, i, i, 0, ADD_VALUES);

    this->PartialAssemble();
}

/**
 * Destructor.
 */
Matrix::~Matrix() {
    this->Destroy();
}

/**
 * Assemble the matrix. This function should be
 * called after all values have been set,
 * and before the matrix is "used".
 */
void Matrix::Assemble() {
    MatAssemblyBegin(this->petsc_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->petsc_mat, MAT_FINAL_ASSEMBLY);
}

/**
 * A partial assemble. This assemble is designed to allow
 * us to switch from ADDING elements to INSERTING elements,
 * and vice versa. This assemble must be followed by a
 * regular 'Assemble()' before more advanced operations
 * (than element insertion) can be conducted.
 */
void Matrix::PartialAssemble() {
    MatAssemblyBegin(this->petsc_mat, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(this->petsc_mat, MAT_FLUSH_ASSEMBLY);
}

/**
 * Checks if this matrix contains any NaN or
 * Inf elements. If arguments are provided, the
 * indices of the first NaN/Inf element are returned.
 *
 * I: Row of first NaN/Inf value (if present).
 * J: Column of first NaN/Inf value (if present).
 */
bool Matrix::ContainsNaNOrInf(len_t *I, len_t *J) {
    const PetscInt *cols;
    const PetscScalar *vals;
    PetscInt ncols;

    for (PetscInt i = 0; i < this->m; i++) {
        MatGetRow(this->petsc_mat, i, &ncols, &cols, &vals);

        for (PetscInt j = 0; j < ncols; j++) {
            if (isnan(vals[j]) || isinf(vals[j])) {
                if (I != nullptr) *I = i;
                if (J != nullptr) *J = cols[j];

                return true;
            }
        }
    }

    return false;
}

/**
 * Destroy this matrix.
 */
void Matrix::Destroy() {
    if (this->allocated)
        MatDestroy(&this->petsc_mat);
}

/**
 * Transform this matrix according to
 *
 *   B = l*A*r,
 *
 * where 'A' is the matrix before the transformation,
 * and 'l' and 'r' are transformation vectors.
 *
 * l: Left vector to transform with (may be 'nullptr').
 * r: Right vector to transform with (may be 'nullptr').
 */
void Matrix::DiagonalScale(Vec l, Vec r) {
    MatDiagonalScale(this->petsc_mat, l, r);
}

/**
 * Returns the index of the first row owned by this
 * CPU as well as the number of rows owned by it.
 *
 * ir: First row of matrix owned by this CPU.
 * nr: Number of rows owned by this CPU.
 *
 * Both parameters contain the desired value on return.
 */
void Matrix::GetOwnershipRange(PetscInt *ir, PetscInt *nr) {
    MatGetOwnershipRange(this->petsc_mat, ir, nr);
    *nr = *nr - *ir;
}

/**
 * Returns the matrix element with the given indices.
 */
real_t Matrix::GetElement(const PetscInt i, const PetscInt j) {
    PetscScalar v;
    MatGetValues(this->petsc_mat, 1, &i, 1, &j, &v);

    return v;
}

/**
 * Returns the requested matrix row in the given vector.
 */
void Matrix::GetRow(const PetscInt i, PetscScalar *v) {
    const PetscInt n = this->n;
    PetscInt *idx = new PetscInt[n];

    for (PetscInt i = 0; i < n; i++)
        idx[i] = i;

    GetRow(i, idx, v);

    delete [] idx;
}
void Matrix::GetRow(const PetscInt i, const PetscInt *j, PetscScalar *v) {
    MatGetValues(this->petsc_mat, 1, &i, n, j, v);
}


/**
 * Returns the requested matrix column in the given vector.
 */
void Matrix::GetColumn(const PetscInt j, PetscScalar *v) {
    PetscInt *idx = new PetscInt[this->m];

    for (PetscInt k = 0; k < this->m; k++)
        idx[k] = k;

    GetColumn(j, idx, v);

    delete [] idx;
}
void Matrix::GetColumn(const PetscInt j, const PetscInt *i, PetscScalar *v) {
    MatGetValues(this->petsc_mat,this->m, i, 1, &j, v);
}



/**
 * Returns the absolute value of the maximum element
 * of each matrix row.
 *
 * v: If provided, puts maxima of each row in this vector.
 *    Must have 'm' elements. If not provided, this function
 *    instead returns a newly allocated vector with 'm'
 *    elements.
 */
real_t *Matrix::GetRowMaxAbs() {
    real_t *v = new real_t[this->m];
    GetRowMaxAbs(v);
    return v;
}
void Matrix::GetRowMaxAbs(real_t *v) {
    Vec s;

    VecCreateSeq(PETSC_COMM_WORLD, this->m, &s);
    VecAssemblyBegin(s);
    VecAssemblyEnd(s);

    MatGetRowMaxAbs(this->petsc_mat, s, NULL);

    PetscInt *idx = new PetscInt[this->m];
    for (PetscInt i = 0; i < this->m; i++)
        idx[i] = i;

    VecGetValues(s, this->m, idx, v);
    VecDestroy(&s);

    delete [] idx;
}

/**
 * Returns the number of non-zero elements in
 * the matrix.
 */
len_t Matrix::GetNNZ() {
    MatInfo info;
    MatGetInfo(this->petsc_mat, MAT_GLOBAL_MAX, &info);

    return info.nz_used;
}

/**
 * Forms the matrix
 * 
 *    A = I - dt * A
 *
 * where I is the unit matrix, A is this matrix and
 * dt is a scalar.
 */
void Matrix::IMinusDtA(const PetscScalar dt) {
    PetscScalar DT = -dt;
    MatScale(this->petsc_mat, DT);
    MatShift(this->petsc_mat, 1.0);
}


/**
 * Multiply this matrix with the the given vector.
 *
 * n: Number of elements in the given vector.
 * f: Vector to multiply with.
 *
 * RETURNS a newly allocated array with n elements
 * containing the result of the multiplication.
 */
real_t *Matrix::Multiply(const len_t n, const real_t *f) {
    if (n != (len_t)this->n)
        throw MatrixException(
            "The given vector has wrong dimensions. Expected vector "
            "to have %llu elements, but it has %llu elements.",
            this->n, n
        );

    Vec f_v, Af_v;

    VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, n, f, &f_v);
    VecCreateSeq(PETSC_COMM_WORLD, this->m, &Af_v);

    VecAssemblyBegin(f_v);  VecAssemblyEnd(f_v);
    VecAssemblyBegin(Af_v); VecAssemblyEnd(Af_v);

    MatMult(this->petsc_mat, f_v, Af_v);

    real_t *Af;
    VecGetArray(Af_v, &Af);

    real_t *r = new real_t[this->m];
    for (PetscInt i = 0; i < this->m; i++)
        r[i] = Af[i];

    VecDestroy(&Af_v);
    VecDestroy(&f_v);

    return r;
}

/**
 * Print basic information about the PETSc matrix.
 */
void Matrix::PrintInfo() {
    MatInfo info;
    MatGetInfo(this->petsc_mat, MAT_GLOBAL_MAX, &info);

    cout << ":: MATRIX INFORMATION" << endl;
    cout << "   Matrix size:            " << this->m << " x " << this->n << endl;
    cout << "   Block size:             " << info.block_size << endl;
    cout << "   Non-zeros" << endl;
    cout << "     allocated:            " << info.nz_allocated << endl;
    cout << "     used:                 " << info.nz_used << endl;
    cout << "     unneeded:             " << info.nz_unneeded << endl;
    cout << "   Memory allocated:       " << info.memory << endl;
    cout << "   Matrix assemblies:      " << info.assemblies << endl;
    cout << "   Number of mallocs:      " << info.mallocs << endl;
    cout << "   Fill ratio for LU/ILU" << endl;
    cout << "     given:                " << info.fill_ratio_given << endl;
    cout << "     needed:               " << info.fill_ratio_needed << endl;
    cout << "   # mallocs during fact.: " << info.factor_mallocs << endl;
}

/**
 * Reset the row and column offsets.
 */
void Matrix::ResetOffset() {
    this->rowOffset = 0;
    this->colOffset = 0;
}

/**
 * Sets one element in the matrix.
 *
 * irow:        Element row.
 * icol:        Element column.
 * v:           Value to write to element.
 * insert_mode: Either 'INSERT_VALUES' or 'ADD_VALUES'.
 */
void Matrix::SetElement(
    const PetscInt irow, const PetscInt icol,
    const PetscScalar v, InsertMode insert_mode
) {
    if(v!=0)
        MatSetValue(this->petsc_mat, this->rowOffset+irow, this->colOffset+icol, v, insert_mode);
}

/**
 * Sets the offset with which elements should be set in
 * the matrix.
 */
void Matrix::SetOffset(const PetscInt rOff, const PetscInt cOff) {
    this->rowOffset = rOff;
    this->colOffset = cOff;
}

/**
 * Sets the values of one row of the matrix.
 */
void Matrix::SetRow(
	PetscInt irow, const PetscInt ncol,
	PetscInt *icol, const PetscScalar *v,
	InsertMode insert_mode
) {
    // Apply offsets
    irow += this->rowOffset;
    for(len_t i=0; i<ncol; i++)
        icol[i] += this->colOffset;
	MatSetValues(this->petsc_mat, 1, &irow, ncol, icol, v, insert_mode);
    // Reset offsets
    for(len_t i=0; i<ncol; i++)
        icol[i] -= this->colOffset;
    irow -= this->rowOffset;
}

/**
 * Print the matrix to STDOUT in the given format.
 *
 * format: Format of output.
 */
void Matrix::View(const enum view_format format, const string& filename) {
    PetscViewer viewer;

    if (format == NON_ZERO_STRUCT)
        viewer = PETSC_VIEWER_DRAW_WORLD;
    else if (format == BINARY_MATLAB) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    } else {
        viewer = PETSC_VIEWER_STDOUT_SELF;
        
        PetscViewerPushFormat(viewer, (PetscViewerFormat)format);
    }

    MatView(this->petsc_mat, viewer);

    if (format == BINARY_MATLAB)
        PetscViewerDestroy(&viewer);
}

/**
 * Zero the non-zero entries. This routine
 * retains the non-zero structure of the matrix, though,
 * and should be called before rebuilding the matrix.
 */
void Matrix::Zero(bool keepNzStructure) {
    if(!keepNzStructure)
        MatSetOption(this->petsc_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE);
    else
        MatSetOption(this->petsc_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    
    MatZeroEntries(this->petsc_mat);
}

/**
 * Zero the non-zero entries of the matrix on the
 * given rows. This routine maintains the non-zero
 * structure of the matrix.
 *
 * n: Number of rows to zero.
 * i: Indices of rows to zero.
 */
void Matrix::ZeroRows(const PetscInt n, const PetscInt i[]) {
	MatZeroRows(this->petsc_mat, n, i, 0.0, nullptr, nullptr);
}

/**
 * Sets diagonal entries of the matrix to a constant value.
 * 
 * n: Number of rows to set to constant
 * i: Indices of rows to set
 * v: The constant value that the diagonal will take 
 */
void Matrix::SetDiagonalConstant(const PetscInt n, const PetscInt i[], const PetscReal v) {
    for(PetscInt it=0; it<n; it++)
        MatSetValue(this->petsc_mat, this->rowOffset+i[it], this->colOffset+i[it], v, INSERT_VALUES);
}
