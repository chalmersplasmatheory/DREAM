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
#include "FVM/QuantityData.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 */
QuantityData::QuantityData(Grid *grid)
    : grid(grid) {

    this->nElements = grid->GetNCells();
    AllocateData();
}

/**
 * Destructor.
 */
QuantityData::~QuantityData() {
    for (auto it = store.begin(); it != store.end(); it++)
        delete [] *it;

    delete [] data;
}

/**
 * Allocate new memory for the temporary data storage
 * of this object. Note that is method does NOT handle
 * cleanup of the previously stored data.
 */
void QuantityData::AllocateData() {
    this->data = new real_t[this->nElements];
    this->idxVec = new PetscInt[this->nElements];

    for (len_t i = 0; i < nElements; i++)
        this->idxVec[i] = (PetscInt)i;
}

/**
 * Save the currently stored data as a new time step.
 * This means that the saved data can be accessed later
 * by requesting the data at the given time step.
 *
 * t: Time to associate with the data.
 */
void QuantityData::SaveStep(const real_t t) {
    times.push_back(t);
    store.push_back(this->data);

    AllocateData();
}

/**
 * Store data from the given PETSc vector into the temporary
 * data store of this 'QuantityData' object.
 *
 * vec:    PETSc vector to copy data from (usually solution vector).
 * offset: Index of first element in the given vector to copy.
 */
void QuantityData::Store(Vec& vec, const len_t offset) {
    if ((len_t)idxVec[0] != offset) {
        for (len_t i = 0; i < nElements; i++)
            idxVec[i] = (PetscInt)(offset + i);
    }

    VecGetValues(vec, (PetscInt)nElements, idxVec, this->data);
}

/**
 * Store data from the given array to into the temporary
 * data store of this 'QuantityData' object.
 *
 * offset: Index of first element in the given vector to copy.
 * vec:    Array to copy data from.
 */
void QuantityData::Store(const real_t *vec, const len_t offset) {
    for (len_t i = 0; i < nElements; i++)
        this->data[i] = vec[offset+i];
}

/**
 * Save the contents of this QuantityData object to a
 * file using the given SFile object.
 *
 * sf:       softlib SFile object to use for saving data.
 * name:     Name of variable in object.
 * path:     Path in SFile to save variable to.
 * saveMeta: If 'true', saves time and coordinate grids along with the variable
 *           data. In this case, 'name' is interpreted as a group name instead
 *           and will contain at least the variables 't' (time) and 'x' (data),
 *           in addition to any coordinate grids (e.g. 'r', 'p', 'xi' etc.).
 */
void QuantityData::SaveSFile(
    SFile *sf, const string& name, const string& path,
    bool saveMeta
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

        // Save time grid
        const real_t *t = this->times.data();
        sf->WriteList(group + "t", t, nt);

        // Write grids
        if (nr > 1) {
            const real_t *r = this->grid->GetRadialGrid()->GetR();
            sf->WriteList(group + "r", r, nr);
        }

        // XXX Here we assume that all momentum grids are the same
        if (np1 > 1) {
            const real_t *p1 = this->grid->GetMomentumGrid(0)->GetP1();
            const string p1name = this->grid->GetMomentumGrid(0)->GetP1Name();

            sf->WriteList(group + p1name, p1, np1);
        }

        if (np2 > 1) {
            const real_t *p2 = this->grid->GetMomentumGrid(0)->GetP2();
            const string p2name = this->grid->GetMomentumGrid(0)->GetP2Name();

            sf->WriteList(group + p2name, p2, np2);
        }
    }

    sfilesize_t ndims = 1;
    sfilesize_t dims[4] = {nt,0,0,0};
    if (nr) dims[ndims++] = nr; 
    if (np2) dims[ndims++] = np2;
    if (np1) dims[ndims++] = np1;

    // Compute number of elements
    len_t nel = dims[0];
    for (len_t i = 1; i < ndims; i++)
        nel *= dims[i];

    // Save data
    if (ndims == 2) {
        // The reason we handle 2D arrays separately is
        // that we can save them without first copying the data
        sf->WriteArray(group + dname, this->store.data(), dims[0], dims[1]);
    } else {    // 1, 3 and 4 dimensions
        // For these arrays, we need to copy the data
        // before saving it...
        real_t *data = new real_t[nel];
        for (len_t i = 0; i < nt; i++) {
            for (len_t j = 0; j < nElements; j++)
                data[i*nElements + j] = this->store[i][j];
        }

        sf->WriteMultiArray(group + dname, data, ndims, dims);

        delete [] data;
    }
}

