/**
 * Implementation of the 'QuantityData' class, which is used to provide an
 * interface for storing selected point in time series of data.
 *
 * In DREAM, this class is used for both temporary and persistent storage of
 * quantities solved for in the 'EquationSystem' object. Each quantity has a
 * 'QuantityData' object belonging to it. Whenever the equation system has
 * been solved, the solution is written to the corresponding 'QuantityData'
 * object, which holds the data temporarily. If the solution in a particular
 * time step should be saved (and, i.e., be part of the final code output),
 * a call to 'SaveStep()' should be made. This automatically records the data
 * and makes it possible to retrieve it later. Otherwise, a call to 'Store()'
 * overwrites the data currently stored (temporarily) in the object. This makes
 * it possible to store only a select set of time steps, thus saving memory.
 */

#include <string>
#include "FVM/FVMException.hpp"
#include "FVM/QuantityData.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 *
 * grid:       Grid on which the quantity is defined.
 * nMultiples: The number of quantities stored within this quantity.
 *             This can be used to store the density of all ions as a
 *             single UnknownQuantity, each defined on a radial grid.
 */
QuantityData::QuantityData(Grid *grid, const len_t nMultiples, enum FVM::fluxGridType fgt)
    : grid(grid), fluxGridType(fgt) {

    len_t n = 1;
    switch (fgt) {
        case FVM::FLUXGRIDTYPE_DISTRIBUTION: n = grid->GetNCells(); break;
        case FVM::FLUXGRIDTYPE_RADIAL: n = grid->GetNCells_fr(); break;
        case FVM::FLUXGRIDTYPE_P1: n = grid->GetNCells_f1(); break;
        case FVM::FLUXGRIDTYPE_P2: n = grid->GetNCells_f2(); break;

        default:
            throw FVM::FVMException("QuantityData: Unrecognized flux grid type specified: %d.", fgt);
    }

    this->nElements  = n * nMultiples;
    this->nMultiples = nMultiples;

    AllocateData();
}

/**
 * Destructor.
 */
QuantityData::~QuantityData() {
    for (auto it = store.begin(); it != store.end(); it++)
        delete [] *it;

    delete [] olddata[0];
    delete [] olddata;
    delete [] data;
    delete [] idxVec;
}

/**
 * Allocate new memory for the temporary data storage
 * of this object.
 */
void QuantityData::AllocateData() {
    this->data = new real_t[this->nElements];

    if (N_SAVE_OLD_STEPS == 0)
        throw FVM::FVMException("QuantityData: The parameter 'N_SAVE_OLD_STEPS' must be at least 1.");

    this->olddata = new real_t*[N_SAVE_OLD_STEPS];
    this->olddata[0] = new real_t[N_SAVE_OLD_STEPS*this->nElements];
    for (len_t i = 1; i < N_SAVE_OLD_STEPS; i++)
        this->olddata[i] = this->olddata[i-1] + this->nElements;

    this->oldtime = new real_t[N_SAVE_OLD_STEPS];
    this->idxVec = new PetscInt[this->nElements];

    this->nOldSaved = 0;

    for (len_t i = 0; i < nElements; i++)
        this->data[i] = 0;
    for (len_t i = 0; i < N_SAVE_OLD_STEPS*nElements; i++)
        this->olddata[0][i] = 0;
    for (len_t i = 0; i < nElements; i++)
        this->idxVec[i] = (PetscInt)i;
}

/**
 * Returns 'true' if we can roll back the saved data by
 * at least one time step.
 */
bool QuantityData::CanRollbackSaveStep() const {
    return (this->nOldSaved > 0);
}

/**
 * Save the currently stored data as a new time step.
 * This means that the saved data can be accessed later
 * by requesting the data at the given time step.
 *
 * t:        Time to associate with the data.
 * trueSave: If 'true', saves the current data to the 'store'
 *           array, in addition to storing it in the 'olddata'
 *           buffer.
 */
void QuantityData::SaveStep(const real_t t, bool trueSave) {
    // Move old data...
    len_t nMin = min(this->nOldSaved, N_SAVE_OLD_STEPS-1);
    for (len_t i = nMin; i > 0; i--) {
        for (len_t j = 0; j < this->nElements; j++)
            this->olddata[i][j] = this->olddata[i-1][j];

        this->oldtime[i] = this->oldtime[i-1];
    }

    if (this->nOldSaved < N_SAVE_OLD_STEPS)
        this->nOldSaved++;

    // Copy previous solution...
    for (len_t i = 0; i < this->nElements; i++)
        this->olddata[0][i] = this->data[i];
    
    this->oldtime[0] = t;

    // Copy to true 'store' array
    if (trueSave) {
        real_t *v = new real_t[this->nElements];
        for (len_t i = 0; i < nElements; i++)
            v[i] = this->olddata[0][i];

        times.push_back(this->oldtime[0]);
        store.push_back(v);

        // If we do a 'true' save, we cannot roll back the solution anymore
        // (this is just an imposed limitation; rolling back stored solutions
        // would require a bit more coding, and as of this writing that
        // functionality is not needed)
        this->nOldSaved = 0;
    }
}

/**
 * Roll back a previously saved time step. This method is
 * the inverse of the method 'SaveStep()' with 'trueSave = false'.
 */
void QuantityData::RollbackSaveStep() {
    if (!CanRollbackSaveStep())
        throw FVM::FVMException("QuantityData: Cannot roll back previous time step.");

    for (len_t i = 0; i < this->nElements; i++)
        this->data[i] = this->olddata[0][i];

    len_t nMin = min(this->nOldSaved, this->N_SAVE_OLD_STEPS);
    for (len_t i = 0; i < nMin-1; i++) {
        for (len_t j = 0; j < this->nElements; j++)
            this->olddata[i][j] = this->olddata[i+1][j];

        this->oldtime[i] = this->oldtime[i+1];
    }

    this->nOldSaved--;
}

/**
 * Store data from the given PETSc vector into the temporary
 * data store of this 'QuantityData' object.
 *
 * vec:           PETSc vector to copy data from (usually solution vector).
 * offset:        Index of first element in the given vector to copy.
 * mayBeConstant: Indicates that the data might have not changed from the
 *                previous iteration and may thus warrant skipping 'Rebuild()'
 *                in certain external objects depending on this data.
 */
void QuantityData::Store(Vec& vec, const len_t offset, bool mayBeConstant) {
    if (mayBeConstant) {
        // Check if the current and given data are
        // exactly equal (i.e. constant)
        Vec tv;
        VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nElements, this->data, &tv);

        PetscBool eq;
        VecEqual(vec, tv, &eq);

        VecDestroy(&tv);

        if (eq == PETSC_TRUE) {
            this->hasChanged = false;
            return;
        }
    }

    this->hasChanged = true;

    if ((len_t)idxVec[0] != offset) {
        for (len_t i = 0; i < nElements; i++)
            idxVec[i] = (PetscInt)(offset + i);
    }

    VecGetValues(vec, (PetscInt)nElements, idxVec, this->data);
}

/**
 * Store data from the given array into the temporary
 * data store of this 'QuantityData' object.
 *
 * offset:        Index of first element in the given vector to copy.
 * vec:           Array to copy data from.
 * mayBeConstant: Indicates that the data might have not changed from the
 *                previous iteration and may thus warrant skipping 'Rebuild()'
 *                in certain external objects depending on this data.
 */
void QuantityData::Store(const real_t *vec, const len_t offset, bool mayBeConstant) {
    if (mayBeConstant) {
        // Check if the current and given data are
        // exactly equal (i.e. constant)
        real_t s = 0;
        for (len_t i = 0; i < nElements; i++)
            s += this->data[i] - vec[offset+i];

        // Equal?
        if (s == 0) {
            this->hasChanged = false;
            return;
        }
    }

    // Data is updated
    this->hasChanged = true;

    for (len_t i = 0; i < nElements; i++)
        this->data[i] = vec[offset+i];
}

/**
 * Store data from the given 2D-array into the temporary
 * data store of this 'QuantityData' object.
 *
 * m, n:          Dimensions of 'vec'.
 * vec:           Data to store.
 * mayBeConstant: If true, indicates that the data may not
 *                have changed since the last iteration.
 *
 * XXX here we assume that all momentum grids have the same
 * sizes at all radii
 */
void QuantityData::Store(
    const len_t m, const len_t n,
    const real_t *const* vec, bool mayBeConstant
) {
    if (mayBeConstant) {
        // Check if the current and given data are
        // exactly equal (i.e. constant)
        real_t s = 0;
        for (len_t i = 0; i < m; i++)
            for (len_t j = 0; j < n; j++)
                s += this->data[i*n + j] - vec[i][j];

        // Equal?
        if (s == 0) {
            this->hasChanged = false;
            return;
        }
    }

    // Data has changed
    this->hasChanged = true;

    for (len_t i = 0; i < m; i++) {
        for (len_t j = 0; j < n; j++)
            this->data[i*n + j] = vec[i][j];
    }
}

/**
 * Store an array of zeros and return a pointer to the
 * data array (this can be used for first making room
 * for the data and _then_ writing to the array).
 */
real_t *QuantityData::StoreEmpty() {
    return this->data;
}

/**
 * Save the contents of this QuantityData object to a
 * file using the given SFile object.
 *
 * sf:          softlib SFile object to use for saving data.
 * name:        Name of variable in object.
 * path:        Path in SFile to save variable to.
 * description: String describing this quantity.
 * saveMeta:    If 'true', saves time and coordinate grids along with the variable
 *              data. In this case, 'name' is interpreted as a group name instead
 *              and will contain at least the variables 't' (time) and 'x' (data),
 *              in addition to any coordinate grids (e.g. 'r', 'p', 'xi' etc.).
 */
void QuantityData::SaveSFile(
    SFile *sf, const string& name, const string& path,
    const string& description, bool saveMeta
) {
    const len_t
        nt  = this->times.size(),
        nr  = this->grid->GetNr(),
        // XXX Here we assume that all momentum grids are the same
        np1 = this->grid->GetMomentumGrid(0)->GetNp1(),
        np2 = this->grid->GetMomentumGrid(0)->GetNp2();

    string group = path + "/";
    string dname = name;

    // Should we save "metadata", i.e. grid data?
    if (saveMeta) {
        group += name + "/";
        dname = "x";

        sf->CreateStruct(group);

        // Save time grid
        const real_t *t = this->times.data();
        sf->WriteList(group + "t", t, nt);

        // Write grids
        if (nr > 1) {
            if (this->fluxGridType == FLUXGRIDTYPE_RADIAL)
                sf->WriteList(group + "r", this->grid->GetRadialGrid()->GetR_f(), nr+1);
            else
                sf->WriteList(group + "r", this->grid->GetRadialGrid()->GetR(), nr);
        }

        // XXX Here we assume that all momentum grids are the same
        if (np1 > 1) {
            const string p1name = this->grid->GetMomentumGrid(0)->GetP1Name();

            if (this->fluxGridType == FLUXGRIDTYPE_P1)
                sf->WriteList(group + p1name, this->grid->GetMomentumGrid(0)->GetP1_f(), np1+1);
            else
                sf->WriteList(group + p1name, this->grid->GetMomentumGrid(0)->GetP1(), np1);
        }

        if (np2 > 1) {
            const string p2name = this->grid->GetMomentumGrid(0)->GetP2Name();

            if (this->fluxGridType == FLUXGRIDTYPE_P2)
                sf->WriteList(group + p2name, this->grid->GetMomentumGrid(0)->GetP2_f(), np2+1);
            else
                sf->WriteList(group + p2name, this->grid->GetMomentumGrid(0)->GetP2(), np2);
        }
    }

    sfilesize_t ndims = 1;
    sfilesize_t dims[5] = {nt,0,0,0,0};

    if (this->nMultiples > 1) dims[ndims++] = this->nMultiples;
    //if (nr > 1 || np2 > 1 || np1 > 1) dims[ndims++] = nr;

    // Always include radial dimension
    if (this->fluxGridType == FLUXGRIDTYPE_RADIAL)
        dims[ndims++] = nr+1;
    else dims[ndims++] = nr;

    if (np2 > 1 || np1 > 1) {
        if (this->fluxGridType == FLUXGRIDTYPE_P2)
            dims[ndims++] = np2+1;
        else dims[ndims++] = np2;

        if (this->fluxGridType == FLUXGRIDTYPE_P1)
            dims[ndims++] = np1+1;
        else dims[ndims++] = np1;
    }

    // Compute number of elements
    len_t nel = dims[0];
    for (len_t i = 1; i < ndims; i++)
        nel *= dims[i];

    // Save data (since it is not stored contiguously in
    // memory, we need to copy it to a new, temporary
    // array first)
    real_t *data = new real_t[nel];
    for (len_t i = 0; i < nt; i++) {
        for (len_t j = 0; j < nElements; j++)
            data[i*nElements + j] = this->store[i][j];
    }

    sf->WriteMultiArray(group + dname, data, ndims, dims);

    if (!description.empty())
        sf->WriteAttribute_string(group + dname, "description", description);

    delete [] data;
}

/**
 * Set the initial value of the specified unknown quantity. If
 * the initial value has previously been specified, it is overwritten.
 * Note that the data is *copied* and, hence, this object does not
 * take over responsibility for freeing the memory occupied by 'val'.
 *
 * val: Initial value of the quantity. If 'nullptr', sets all
 *      elements to zero in the initial value.
 * t0:  Initial time.
 */
void QuantityData::SetInitialValue(const real_t *val, const real_t t0) {
    bool hasValue = this->HasInitialValue();

    const real_t *init = val;
    if (val == nullptr) {
        real_t *_init = new real_t[nElements];
        for (len_t i = 0; i < nElements; i++)
            _init[i] = 0;

        init = _init;
    }

    this->Store(init);
    
    // Save to 'store' array (which is written to final output)
    if (!hasValue)
        this->SaveStep(t0, true);
    else {
        this->SaveStep(t0, false);

        // Reset nOldSaved since we didn't push new data just now,
        // but rather just replace it.
        this->nOldSaved = 0;

        // Overwrite previously stored value
        real_t *iv = this->store[0];
        this->times[0] = t0;

        for (len_t i = 0; i < nElements; i++)
            iv[i] = init[i];
    }

    if (val == nullptr)
        delete [] init;
}

