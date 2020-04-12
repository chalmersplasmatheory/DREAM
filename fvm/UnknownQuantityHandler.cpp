/**
 * Implementation of an unknown quantity handler (a collection
 * of unknowns).
 */

#include <string>
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
 * Returns the most recent data for the specified
 * unknown quantity.
 *
 * qty: ID of quantity to get data of.
 */
real_t *UnknownQuantityHandler::GetUnknownData(const len_t qty) {
    return unknowns[qty]->GetData();
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
 * Add an unknown quantity to the equation system.
 *
 * name: Name of unknown quantity.
 * grid: Grid on which the quantity is defined.
 */
len_t UnknownQuantityHandler::InsertUnknown(const string& name, FVM::Grid *grid) {
    unknowns.push_back(new UnknownQuantity(name, grid));

    // Return ID of quantity
    return (unknowns.size()-1);
}

/**
 * Store the given vector to the specified list of unknowns.
 *
 * unk: List of unknown IDs to store elements of vector to.
 * v:   Vector containing the data to store.
 */
void UnknownQuantityHandler::Store(vector<len_t> &unk, Vec &v) {
    len_t offset = 0;
    for (auto it = unk.begin(); it != unk.end(); it++) {
        unknowns[*it]->Store(v, offset);
        offset += unknowns[*it]->NumberOfElements();
    }
}

