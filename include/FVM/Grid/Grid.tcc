/**
 * Implementation of the compound 'Grid' object,
 * which represents a full computational grid.
 */

#include "FVM/Grid/Grid1D.hpp"


/**
 * Generate this Grid from a number of 1D
 * Grid dimensions.
 */
template<int N, typename ... Args>
TQS::FVM::Grid<N>(Args&& ... args) {
    static_assert(sizeof(Args)==N, "Invalid number of Grid1D objects given to the 'Grid' object.");

    this->insert_dimensions(args...);
}

/**
 * Copy constructor.
 */
template<int N>
TQS::FVM::Grid<N>(Grid *g) {
}


/************************************
 * PRIVATE METHODS                  *
 ************************************/
/**
 * Insert a grid dimension into this grid.
 *
 * g: Grid dimension object to insert.
 */
template<int N>
void TQS::FVM::Grid<N>::insert_dimensions(Grid1D *g) {
    this->dimensions.push_back(g);
}
template<int N, typename ... Args>
void TQS::FVM::Grid<N>::insert_dimensions(Grid1D *g, Args&& ... args) {
    this->dimensions.push_back(g);
    this->insert_dimensions(args...);
}

/**
 * Converts a list of indices to a linear index.
 */
len_t get_index(const len_t i) const {
    return i;
}
template<int N, typename ... Args>
len_t TQS::FVM::Grid<N>::get_index(Args&& ... args, const len_t i) const {
    len_t nx = this->dimensions[N-sizeof(Args)+1].size();
    return ((get_index(args) * nx) + i);
}


/************************************
 * PUBLIC METHODS                   *
 ************************************/
/**
 * Rebuild this grid for the specified time point.
 * Note that this does NOT force a grid rebuild,
 * and so does NOT rebuild static (i.e. time-independent)
 * grids. If at least one of the Grid dimensions
 * needs to be rebuilt, all quantities owned by
 * this grid object will also be rebuilt (i.e.
 * cell volumes etc.)
 *
 * t: Time for which to rebuild the grid.
 *
 * RETURNS true if the grid was rebuilt; false otherwise.
 */
bool TQS::FVM::Grid<N>::Rebuild(const real_t t) {
    bool rebuildNeeded = false;
    for (int i = 0; i < N; i++) {
        rebuildNeeded &= this->dimensions[i].Rebuild(t);
    }
}

/**
 * Returns the volume of the given cell. This function
 * takes a list of indices into the grid as argument.
 */
template<int N, typename ... Args>
real_t TQS::FVM::Grid<N>::Volume(Args&& ... args) const {
    static_assert(sizeof(Args)<=N, "Too many indices given.");
    static_assert(sizeof(Args)==N, "Not enough indices given.");

    len_t i = get_index(args...);
    return this->volumes[i];
}

