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

#include "DREAM/QuantityData.hpp"


using namespace DREAM;


/**
 * Constructor.
 */
QuantityData::QuantityData(FVM::Grid *grid)
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
 * offset: Index of first element in the given vector to copy.
 * vec:    PETSc vector to copy data from (usually solution vector).
 */
void QuantityData::Store(const len_t offset, Vec& vec) {
    if (idxVec[0] != offset) {
        for (len_t i = 0; i < nElements; i++)
            idxVec[i] = (PetscInt)(offset + i);
    }

    VecGetValues(vec, (PetscInt)nElements, idxVec, this->data);
}

