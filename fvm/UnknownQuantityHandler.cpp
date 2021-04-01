/**
 * Implementation of an unknown quantity handler (a collection
 * of unknowns).
 */

#include <string>
#include <softlib/SFile.h>
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM::FVM;
using namespace std;

/**
 * Constructor.
 */
UnknownQuantityHandler::UnknownQuantityHandler() { }

/**
 * Destructor.
 */
UnknownQuantityHandler::~UnknownQuantityHandler() {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++)
        delete (*it);
}

/**
 * Returns the unknown data in a single, contiguously stored vector for
 * all the specified unknowns.
 *
 * nontrivial_unknowns: List of unknowns to get data for
 *                      (these are usually the "non-trivial" unknowns that
 *                      appear in the equation system to solve.
 */
const real_t *UnknownQuantityHandler::GetLongVector(const vector<len_t>& nontrivial_unknowns, real_t *vec) {
    return GetLongVector(nontrivial_unknowns.size(), nontrivial_unknowns.data(), vec);
}
const real_t *UnknownQuantityHandler::GetLongVector(const len_t n, const len_t *iuqn, real_t *vec) {
    const len_t size = GetLongVectorSize(n, iuqn);

    if (vec == nullptr)
        vec = new real_t[size];

    for (len_t i = 0, offset = 0; i < n; i++) {
        UnknownQuantity *uqn = unknowns[iuqn[i]];
        const len_t N = uqn->NumberOfElements();
        const real_t *data = uqn->GetData();

        for (len_t j = 0; j < N; j++)
            vec[offset + j] = data[j];

        offset += N;
    }

    return vec;
}

/**
 * Returns the unknown data for the previous time step in a single,
 * contiguously stored vector for all the specified unknowns.
 *
 * nontrivial_unknowns: List of unknowns to get data for.
 *                      (these are usually the "non-trivial" unknowns that
 *                      appear in the equation system solved)
 */
const real_t *UnknownQuantityHandler::GetLongVectorPrevious(const vector<len_t>& nontrivial_unknowns, real_t *vec) {
    return GetLongVectorPrevious(nontrivial_unknowns.size(), nontrivial_unknowns.data(), vec);
}
const real_t *UnknownQuantityHandler::GetLongVectorPrevious(const len_t n, const len_t *iuqn, real_t *vec) {
    const len_t size = GetLongVectorSize(n, iuqn);

    if (vec == nullptr)
        vec = new real_t[size];

    for (len_t i = 0, offset = 0; i < n; i++) {
        UnknownQuantity *uqn = unknowns[iuqn[i]];
        const len_t N = uqn->NumberOfElements();
        const real_t *data = uqn->GetDataPrevious();

        for (len_t j = 0; j < N; j++)
            vec[offset + j] = data[j];

        offset += N;
    }

    return vec;
}

/**
 * Returns the data for all unknowns in a single long vector.
 *
 * vec: Vector to store data in. Must be of the size reported by
 *      'GetLongVectorSizeAll()'. If 'nullptr', new memory is allocated
 *      by this method.
 */
const real_t *UnknownQuantityHandler::GetLongVectorAll(real_t *vec) {
    const len_t size = GetLongVectorSizeAll();

    if (vec == nullptr)
        vec = new real_t[size];

    const len_t n = this->unknowns.size();
    for (len_t i = 0, offset = 0; i < n; i++) {
        UnknownQuantity *uqn = unknowns[i];
        const len_t N = uqn->NumberOfElements();
        const real_t *data = uqn->GetData();

        for (len_t j = 0; j < N; j++)
            vec[offset + j] = data[j];

        offset += N;
    }

    return vec;
}

/**
 * Returns the size of the long vector which would be returned by
 * 'GetLongVector()'. Put differently: return the combined number of
 * elements in the unknowns specified to this routine.
 */
const len_t UnknownQuantityHandler::GetLongVectorSize(const vector<len_t>& nontrivial_unknowns) {
    return GetLongVectorSize(nontrivial_unknowns.size(), nontrivial_unknowns.data());
}
const len_t UnknownQuantityHandler::GetLongVectorSize(const len_t n, const len_t *iuqn) {
    len_t size = 0;
    for (len_t i = 0; i < n; i++)
        size += unknowns[iuqn[i]]->NumberOfElements();

    return size;
}

/**
 * Returns the size of the long vector which contains _all_ unknown
 * quantities.
 */
const len_t UnknownQuantityHandler::GetLongVectorSizeAll() {
    len_t size = 0;
    for (UnknownQuantity *uqn : unknowns)
        size += uqn->NumberOfElements();
    
    return size;
}

/**
 * Returns the length of the time step most recently taken
 * and stored (but not necessarily saved) for the specified
 * unknown quantity.
 */
real_t UnknownQuantityHandler::GetUnknownDataPreviousTime(const len_t qty) {
    return unknowns[qty]->GetPreviousTime();
}

/**
 * Returns the most recent data for the specified
 * unknown quantity.
 *
 * qty: ID of quantity to get data of.
 */
real_t *UnknownQuantityHandler::GetUnknownData(const len_t qty) {
    return unknowns[qty]->GetData();
}
/**
 * name: Name of unknown quantity.
 */
real_t *UnknownQuantityHandler::GetUnknownData(const std::string& name) {
    return GetUnknownData(GetUnknownID(name));
}

/**
 * Returns the data for the specified unknown in the
 * previous time step. Note that 'GetUnknownData()' and
 * 'GetUnknownDataPrevious()' return the same data, corresponding
 * to the most recently stored time step, if no iterations
 * have yet been made for the current time step.
 *
 * qty: ID of quantity to get data of.
 */
real_t *UnknownQuantityHandler::GetUnknownDataPrevious(const len_t qty) {
    return unknowns[qty]->GetDataPrevious();
}
/**
 * name: Name of unknown quantity.
 */
real_t *UnknownQuantityHandler::GetUnknownDataPrevious(const std::string& name) {
    return GetUnknownDataPrevious(GetUnknownID(name));
}

/**
 * Returns the initial data for the specified unknown quantity.
 *
 * qty: ID of quantity to get data of.
 */
real_t *UnknownQuantityHandler::GetUnknownInitialData(const len_t qty) {
    return unknowns[qty]->GetInitialData();
}

/**
 * Returns the ID of the named unknown.
 *
 * name: Name of unknown quantity to get ID of.
 */
len_t UnknownQuantityHandler::GetUnknownID(const string& name) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++) {
        if ((*it)->GetName() == name)
            return (it-unknowns.begin());
    }

    throw FVMException(
        "No unknown quantity with name '%s' exists in the equation system.",
        name.c_str()
    );
}

/**
 * Checks whether an unknown quantity with the given name
 * exists in this UnknownQuantityHandler.
 *
 * name: Name of unknown to look for.
 */
bool UnknownQuantityHandler::HasUnknown(const string& name) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++) {
        if ((*it)->GetName() == name)
            return true;
    }

    return false;
}

/**
 * Add an unknown quantity to the equation system.
 *
 * name: Name of unknown quantity.
 * grid: Grid on which the quantity is defined.
 */
len_t UnknownQuantityHandler::InsertUnknown(const string& name, const string& desc, FVM::Grid *grid, const len_t nMultiples) {
    unknowns.push_back(new UnknownQuantity(name, desc, grid, nMultiples));

    // Return ID of quantity
    return (unknowns.size()-1);
}

/**
 * Store the given vector to the specified list of unknowns.
 *
 * unk:           List of unknown IDs to store elements of vector to.
 * v:             Vector containing the data to store.
 * mayBeConstant: Indicates that the data might have not changed from the
 *                previous iteration and may thus warrant skipping 'Rebuild()'
 *                in certain external objects depending on this data.
 */
void UnknownQuantityHandler::Store(vector<len_t> &unk, Vec &v, bool mayBeConstant) {
    len_t offset = 0;
    for (auto it = unk.begin(); it != unk.end(); it++) {
        unknowns[*it]->Store(v, offset, mayBeConstant);
        offset += unknowns[*it]->NumberOfElements();
    }
}
void UnknownQuantityHandler::Store(vector<len_t> &unk, const real_t *v, bool mayBeConstant) {
    len_t offset = 0;
    for (auto it = unk.begin(); it != unk.end(); it++) {
        unknowns[*it]->Store(v, offset, mayBeConstant);
        offset += unknowns[*it]->NumberOfElements();
    }
}

/**
 * Sets all unknown quantities from the given long vector.
 *
 * vec: Long vector containing data for _all_ unknown quantities.
 */
void UnknownQuantityHandler::SetFromLongVectorAll(const real_t *vec, bool mayBeConstant) {
    len_t offset = 0;
    for (UnknownQuantity *uqn : unknowns) {
        uqn->Store(vec, offset, mayBeConstant);
        offset += uqn->NumberOfElements();
    }
}

/**
 * Rolls back the solution by one time step. This method can
 * only be called once after each 'SaveStep()' has been called.
 */
void UnknownQuantityHandler::RollbackSaveStep() {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++)
        (*it)->RollbackSaveStep();
}

/**
 * Save data for the current time step.
 *
 * t: Time corresponding to the step to save.
 */
void UnknownQuantityHandler::SaveStep(const real_t t, bool trueSave) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++)
        (*it)->SaveStep(t, trueSave);
}

/**
 * Save this list of unknonws to a new file with the given name.
 *
 * filename: Name of file to save data to.
 * saveMeta: If 'true', stores time and coordinate grids along
 *           with the data.
 */
void UnknownQuantityHandler::SaveSFile(
    const string& filename, bool saveMeta
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    this->SaveSFile(sf, "", saveMeta);

    sf->Close();
}

/**
 * Save this list of unknowns to a file using the given SFile object.
 *
 * sf:       SFile object to use for saving data.
 * path:     Path in file to save data to.
 * saveMeta: If 'true', stores time and coordinate grids along
 *           with the data.
 */
void UnknownQuantityHandler::SaveSFile(
    SFile *sf, const string& path, bool saveMeta
) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++)
        (*it)->SaveSFile(sf, path, saveMeta);
}

/**
 * Save the most recent data for this list of unknowns to a
 * new file with the given name.
 *
 * filename: Name of file to save data to.
 * saveMeta: If 'true', stores time and coordinate grids along
 *           with the data.
 */
void UnknownQuantityHandler::SaveSFileCurrent(
    const string& filename, bool saveMeta
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    this->SaveSFileCurrent(sf, "", saveMeta);

    sf->Close();
}

/**
 * Save the most recent data for this list of unknowns to a
 * file using the given SFile object.
 *
 * sf:       SFile object to use for saving data.
 * path:     Path in file to save data to.
 * saveMeta: If 'true', stores time and coordinate grids along
 *           with the data.
 */
void UnknownQuantityHandler::SaveSFileCurrent(
    SFile *sf, const string& path, bool saveMeta
) {
    for (auto it = unknowns.begin(); it != unknowns.end(); it++)
        (*it)->SaveSFileCurrent(sf, path, saveMeta);
}

/**
 * Set the initial value of the specified unknown quantity. If
 * the initial value has previously been specified, it is overwritten.
 *
 * id:  ID of unknown quantity.
 * val: Initial value of the quantity.
 * t0:  Initial time.
 */
void UnknownQuantityHandler::SetInitialValue(
    const string& name, const real_t *val, const real_t t0
) {
    this->SetInitialValue(GetUnknownID(name), val, t0);
}
void UnknownQuantityHandler::SetInitialValue(
    const len_t id, const real_t *val, const real_t t0
) {
    this->unknowns[id]->SetInitialValue(val, t0);
}

