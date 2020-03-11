/**
 * Implementation of the compound 'Grid' object,
 * which represents a full computational grid.
 */

#include "FVM/Grid/GridDimension.hpp"


/**
 * Constructor.
 */
template<len_t N>
TQS::FVM::Grid<N>(const len_t ndims[N]) {
    for (len_t i = 0; i < N; i++)
        this->ndims[i] = ndims[i];
}

/************************************
 * PUBLIC METHODS                   *
 ************************************/
/**
 * Returns the value of the n'th coordinate in the
 * specified cell grid point.
 *
 * n:   Coordinate index.
 * ...: Index tuple specifying the grid point in
 *      which to evaluate the coordinate.
 */
template<len_t N, typename ... Args>
real_t TQS::FVM::Grid<N>::GetX(const len_t n, Args&& ... args) {
    static_assert(sizeof(Args)<=N, "Too many indices given.");
    static_assert(sizeof(Args)==N, "Not enough indices given.");

    return this->get_xn(n, args...);
}

/**
 * Returns the value of the n'th coordinate in the
 * specified flux grid point.
 *
 * n:   Coordinate index.
 * ...: Index tuple specifying the grid point in
 *      which to evaluate the coordinate.
 */
template<len_t N, typename ... Args>
real_t TQS::FVM::Grid<N>::GetX_f(const len_t n, Args&& ... args) {
    static_assert(sizeof(Args)<=N, "Too many indices given.");
    static_assert(sizeof(Args)==N, "Not enough indices given.");

    return this->get_xn(n, args...);
}

/**
 * Rebuild this grid for the specified time point.
 * Note that this does NOT force a grid rebuild,
 * and so does NOT rebuild static (i.e. time-independent)
 * grids.
 *
 * t: Time for which to rebuild the grid.
 *
 * RETURNS true if the grid was rebuilt; false otherwise.
 */
bool TQS::FVM::Grid<N>::Rebuild(const real_t t) {
    return this->rebuild(this->volumes);
}

/**
 * Returns the volume of the specified cell. This function
 * takes a list of indices into the grid as argument.
 */
template<len_t N, typename ... Args>
real_t TQS::FVM::Grid<N>::Volume(Args&& ... args) const {
    static_assert(sizeof(Args)<=N, "Too many indices given.");
    static_assert(sizeof(Args)==N, "Not enough indices given.");

    len_t i = get_index(args...);
    return this->volumes[i];
}

