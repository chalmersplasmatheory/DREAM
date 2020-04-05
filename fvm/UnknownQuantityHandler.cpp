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
        "No unknown quantity with name '%s' exists in the equation system."
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

