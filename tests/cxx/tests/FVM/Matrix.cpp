/**
 * Tests for the 'Matrix' class in the FVM library, focusing on the
 * element insertion routines (buffered 'SetElement()' and 'SetColumn()').
 *
 * All values used in these tests are exactly representable in floating
 * point (multiples of 1/4), so that matrices built via different insertion
 * paths can be compared exactly, independently of summation order.
 */

#include <algorithm>
#include <random>
#include <vector>
#include <petscmat.h>
#include "FVM/Matrix.hpp"
#include "Matrix.hpp"


using namespace DREAMTESTS::FVM;
using namespace std;


// Size of the full (square) matrix
static const PetscInt N = 16;
// Offset and size of the sub-block that entries are written to
static const PetscInt ROW_OFFSET = 3, COL_OFFSET = 5;
static const PetscInt NROWS = 10, NCOLS = 9;
// Number of elements in the long row (exceeding the reserve
// several times over, to exercise buffer growth)
static const PetscInt NLONG = 40;


/**
 * Construct a new matrix of size N-by-N, preallocated as dense
 * (so that no insertion pattern can trigger a preallocation error).
 */
DREAM::FVM::Matrix *Matrix::NewMatrix() {
    return new DREAM::FVM::Matrix(N, N, N);
}

/**
 * Compare the two given (assembled) matrices exactly, both with
 * respect to the non-zero structure and the values.
 *
 * A, B: Matrices to compare.
 * name: Name of the test (for error messages).
 */
bool Matrix::CompareMatrices(
    DREAM::FVM::Matrix *A, DREAM::FVM::Matrix *B, const string& name
) {
    bool success = true;

    for (PetscInt i = 0; i < N && success; i++) {
        PetscInt nA, nB;
        const PetscInt *cA, *cB;
        const PetscScalar *vA, *vB;

        MatGetRow(A->mat(), i, &nA, &cA, &vA);
        MatGetRow(B->mat(), i, &nB, &cB, &vB);

        if (nA != nB) {
            this->PrintError(
                "%s: Number of structural non-zeros differ in row %" PetscInt_FMT
                ": %" PetscInt_FMT " != %" PetscInt_FMT ".",
                name.c_str(), i, nA, nB
            );
            success = false;
        } else {
            for (PetscInt k = 0; k < nA; k++) {
                if (cA[k] != cB[k]) {
                    this->PrintError(
                        "%s: Non-zero structure differs in row %" PetscInt_FMT
                        ": column %" PetscInt_FMT " != %" PetscInt_FMT ".",
                        name.c_str(), i, cA[k], cB[k]
                    );
                    success = false;
                    break;
                } else if (vA[k] != vB[k]) {
                    this->PrintError(
                        "%s: Element (%" PetscInt_FMT ", %" PetscInt_FMT
                        ") differs: %.4f != %.4f.",
                        name.c_str(), i, cA[k],
                        (double)PetscRealPart(vA[k]), (double)PetscRealPart(vB[k])
                    );
                    success = false;
                    break;
                }
            }
        }

        MatRestoreRow(B->mat(), i, &nB, &cB, &vB);
        MatRestoreRow(A->mat(), i, &nA, &cA, &vA);
    }

    return success;
}

/**
 * Check whether element (i,j) (global indices) is part of the
 * non-zero structure of the given (assembled) matrix.
 */
bool Matrix::HasStructuralEntry(
    DREAM::FVM::Matrix *A, const PetscInt i, const PetscInt j
) {
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    bool found = false;

    MatGetRow(A->mat(), i, &ncols, &cols, &vals);
    for (PetscInt k = 0; k < ncols && !found; k++)
        found = (cols[k] == j);
    MatRestoreRow(A->mat(), i, &ncols, &cols, &vals);

    return found;
}

/**
 * Generate a (deterministic) list of matrix entries, given in local
 * (sub-block) coordinates. The list contains exact zeros, rows without
 * any entries, and (unless 'unique' is true) repeated (row,col) pairs,
 * one row with many entries, and entries revisiting earlier rows after
 * other rows have been visited.
 *
 * unique: If true, each (row,col) pair appears at most once in the
 *         list and rows appear in contiguous groups (suitable for
 *         INSERT_VALUES).
 */
vector<struct Matrix::entry> Matrix::GenerateEntries(bool unique) {
    vector<struct entry> entries;
    mt19937 gen(1729);
    uniform_int_distribution<int> dval(-4, 4);
    uniform_int_distribution<PetscInt> dcol(0, NCOLS-1);
    uniform_int_distribution<int> dcount(1, NCOLS);

    auto randval = [&gen,&dval]() { return PetscScalar(dval(gen)) / 4.0; };

    vector<PetscInt> cols(NCOLS);
    for (PetscInt j = 0; j < NCOLS; j++)
        cols[j] = j;

    for (PetscInt i = 0; i < NROWS; i++) {
        // Leave some rows empty
        if (i % 4 == 2)
            continue;

        // Long row with many (repeated) entries
        if (!unique && i == 5) {
            for (PetscInt k = 0; k < NLONG; k++)
                entries.push_back({i, dcol(gen), randval()});
            continue;
        }

        const int n = dcount(gen);
        if (unique) {
            shuffle(cols.begin(), cols.end(), gen);
            for (int k = 0; k < n; k++)
                entries.push_back({i, cols[k], randval()});
        } else {
            for (int k = 0; k < n; k++)
                entries.push_back({i, dcol(gen), randval()});
        }

        // Explicit zero (should never create a structural non-zero)
        if (!unique)
            entries.push_back({i, dcol(gen), 0.0});
    }

    // Revisit earlier rows (only meaningful with ADD_VALUES)
    if (!unique) {
        const PetscInt revisit[] = {0, 7, 1, 1, 5, 0};
        for (PetscInt i : revisit)
            for (PetscInt k = 0; k < 3; k++)
                entries.push_back({i, dcol(gen), randval()});
    }

    return entries;
}

/**
 * Verify that buffered insertion with ADD_VALUES gives exactly the
 * same matrix as unbuffered insertion. Also verifies that unbuffered
 * insertion works correctly after 'EndBufferedSet()', and that a
 * second buffered region (re-using the internal buffers) is correct.
 */
bool Matrix::CheckBufferedAdd() {
    vector<struct entry> entries = GenerateEntries(false);
    DREAM::FVM::Matrix *ref = NewMatrix();
    DREAM::FVM::Matrix *buf = NewMatrix();

    ref->SetOffset(ROW_OFFSET, COL_OFFSET);
    buf->SetOffset(ROW_OFFSET, COL_OFFSET);

    // Region 1: buffered, with a reserve smaller than the longest row
    for (auto &e : entries)
        ref->SetElement(e.row, e.col, e.val, ADD_VALUES);

    buf->BeginBufferedSet(ADD_VALUES, 2);
    for (auto &e : entries)
        buf->SetElement(e.row, e.col, e.val);
    buf->EndBufferedSet();

    // Unbuffered insertion after ending buffered region
    ref->SetElement(4, 3, 1.25, ADD_VALUES);
    buf->SetElement(4, 3, 1.25, ADD_VALUES);

    // Region 2: buffered again, re-using buffers
    buf->BeginBufferedSet(ADD_VALUES);
    for (PetscInt k = 0; k < 10; k++) {
        const PetscInt i = (3*k) % NROWS, j = (5*k) % NCOLS;
        const PetscScalar v = 0.25*(k+1);
        ref->SetElement(i, j, v, ADD_VALUES);
        buf->SetElement(i, j, v);
    }
    buf->EndBufferedSet();

    ref->Assemble();
    buf->Assemble();

    bool success = CompareMatrices(ref, buf, "Buffered ADD_VALUES");

    delete buf;
    delete ref;

    return success;
}

/**
 * Verify that buffered insertion with INSERT_VALUES gives exactly the
 * same matrix as unbuffered insertion.
 */
bool Matrix::CheckBufferedInsert() {
    vector<struct entry> entries = GenerateEntries(true);
    DREAM::FVM::Matrix *ref = NewMatrix();
    DREAM::FVM::Matrix *buf = NewMatrix();

    ref->SetOffset(ROW_OFFSET, COL_OFFSET);
    buf->SetOffset(ROW_OFFSET, COL_OFFSET);

    for (auto &e : entries)
        ref->SetElement(e.row, e.col, e.val, INSERT_VALUES);

    buf->BeginBufferedSet(INSERT_VALUES);
    for (auto &e : entries)
        buf->SetElement(e.row, e.col, e.val, INSERT_VALUES);
    buf->EndBufferedSet();

    ref->Assemble();
    buf->Assemble();

    bool success = CompareMatrices(ref, buf, "Buffered INSERT_VALUES");

    delete buf;
    delete ref;

    return success;
}

/**
 * Verify that entries in the last (and here, only) buffered row are
 * flushed by 'EndBufferedSet()', even though the row never changes.
 */
bool Matrix::CheckBufferedFinalRowFlush() {
    bool success = true;
    DREAM::FVM::Matrix *buf = NewMatrix();
    buf->SetOffset(ROW_OFFSET, COL_OFFSET);

    const PetscInt irow = 2;
    const PetscInt icol[] = {0, 4, 8, 4};
    const PetscScalar v[] = {1.0, -0.5, 2.25, 0.75};

    buf->BeginBufferedSet(ADD_VALUES);
    for (PetscInt k = 0; k < 4; k++)
        buf->SetElement(irow, icol[k], v[k]);
    buf->EndBufferedSet();
    buf->Assemble();

    const PetscInt I = ROW_OFFSET + irow;
    const PetscInt J[] = {COL_OFFSET+0, COL_OFFSET+4, COL_OFFSET+8};
    const PetscScalar expected[] = {1.0, 0.25, 2.25};

    for (PetscInt k = 0; k < 3; k++) {
        const real_t val = buf->GetElement(I, J[k]);
        if (val != expected[k]) {
            this->PrintError(
                "Final buffered row not flushed correctly: element ("
                "%" PetscInt_FMT ", %" PetscInt_FMT ") = %.4f, expected %.4f.",
                I, J[k], (double)val, (double)PetscRealPart(expected[k])
            );
            success = false;
        }
    }

    delete buf;

    return success;
}

/**
 * Verify that 'SetColumn()' inserts the expected values, applying
 * row/column offsets and skipping exact zeros (i.e. not creating
 * structural non-zeros for them), for both ADD_VALUES and INSERT_VALUES.
 */
bool Matrix::CheckSetColumn() {
    bool success = true;

    const PetscInt firstRow = 1, icol = 2, m = 7;
    const PetscScalar vals1[m] = {1.5, 0.0, -0.25, 0.0, 3.0, 0.0, 0.5};
    const PetscScalar vals2[m] = {0.25, 0.0, 0.75, 0.0, 0.0, 0.0, -1.0};

    for (InsertMode mode : {ADD_VALUES, INSERT_VALUES}) {
        const string name = (mode == ADD_VALUES ? "SetColumn ADD_VALUES" : "SetColumn INSERT_VALUES");

        DREAM::FVM::Matrix *ref = NewMatrix();
        DREAM::FVM::Matrix *col = NewMatrix();

        ref->SetOffset(ROW_OFFSET, COL_OFFSET);
        col->SetOffset(ROW_OFFSET, COL_OFFSET);

        for (PetscInt r = 0; r < m; r++)
            ref->SetElement(firstRow+r, icol, vals1[r], mode);
        col->SetColumn(firstRow, icol, m, vals1, mode);

        // Second call on the same column (accumulates/overwrites),
        // and a call with no rows (should be a no-op)
        for (PetscInt r = 0; r < m; r++)
            ref->SetElement(firstRow+r, icol, vals2[r], mode);
        col->SetColumn(firstRow, icol, m, vals2, mode);
        col->SetColumn(firstRow, icol+1, 0, vals1, mode);

        ref->Assemble();
        col->Assemble();

        if (!CompareMatrices(ref, col, name))
            success = false;

        // Explicitly check values, and that zeros did not create
        // structural non-zeros (diagonal elements are always
        // structurally present, so they are skipped)
        const PetscInt J = COL_OFFSET + icol;
        for (PetscInt r = 0; r < m; r++) {
            const PetscInt I = ROW_OFFSET + firstRow + r;
            PetscScalar expected;
            if (mode == ADD_VALUES)
                expected = vals1[r] + vals2[r];
            else
                expected = (vals2[r] != 0 ? vals2[r] : vals1[r]);

            const real_t v = col->GetElement(I, J);
            if (v != expected) {
                this->PrintError(
                    "%s: element (%" PetscInt_FMT ", %" PetscInt_FMT ") = %.4f, expected %.4f.",
                    name.c_str(), I, J, (double)v, (double)PetscRealPart(expected)
                );
                success = false;
            }

            if (vals1[r] == 0 && vals2[r] == 0 && I != J && HasStructuralEntry(col, I, J)) {
                this->PrintError(
                    "%s: zero value created structural non-zero at ("
                    "%" PetscInt_FMT ", %" PetscInt_FMT ").",
                    name.c_str(), I, J
                );
                success = false;
            }
        }

        delete col;
        delete ref;
    }

    return success;
}

/**
 * Run all Matrix tests.
 * Returns 'true' if all tests passed. 'false' otherwise.
 */
bool Matrix::Run(bool) {
    bool success = true;

    if (CheckBufferedAdd())
        this->PrintOK("Buffered SetElement() with ADD_VALUES matches unbuffered insertion.");
    else {
        success = false;
        this->PrintError("Buffered SetElement() with ADD_VALUES does not match unbuffered insertion.");
    }

    if (CheckBufferedInsert())
        this->PrintOK("Buffered SetElement() with INSERT_VALUES matches unbuffered insertion.");
    else {
        success = false;
        this->PrintError("Buffered SetElement() with INSERT_VALUES does not match unbuffered insertion.");
    }

    if (CheckBufferedFinalRowFlush())
        this->PrintOK("EndBufferedSet() flushes the final buffered row.");
    else {
        success = false;
        this->PrintError("EndBufferedSet() does not flush the final buffered row correctly.");
    }

    if (CheckSetColumn())
        this->PrintOK("SetColumn() inserts the expected values.");
    else {
        success = false;
        this->PrintError("SetColumn() does not insert the expected values.");
    }

    return success;
}

