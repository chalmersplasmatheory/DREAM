/**
 * Implementation of the general 'MomentumGrid' class.
 *
 * This grid consists of two arbitrary momentum coordinates, referred to as 'p1'
 * and 'p2'. For each coordinate, this grid contains corresponding cell and flux
 * coordinate grids, as well as the associated cell/flux coordinate step lengths
 * and cell areas.
 *
 * The grid must be generated with an appropriate 'MomentumGridGenerator', which
 * sets the details of the coordinate system.
 *
 * NOTE that this class is abstract, and, although it implements most of the
 * functionality needed by a momentum grid, cannot be used directly. Instead,
 * you must use one of the more specific classes derived from this one.
 */


#include "FVM/Grid/MomentumGrid.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * g:  The generator to use for building this momentum grid.
 * ri: Index of the radial grid point to which this momentum
 *     grid corresponds.
 * rg: Radial grid on which this momentum grid lives.
 * t0: Time at which to initialize the grid.
 */
MomentumGrid::MomentumGrid(MomentumGridGenerator *g, const len_t /*ri*/, const RadialGrid */*rGrid*/, const real_t /*t0*/)
    : generator(g) {

//    g->Rebuild(t0, ri, this, rGrid);
}

/**
 * Destructor.
 */
MomentumGrid::~MomentumGrid() {
    DeallocateP1();
    DeallocateP2();
    DeallocatePAndXi0();

    delete this->generator;
}

/**
 * Deallocation routines.
 */
void MomentumGrid::DeallocateP1() {
    if (np1 == 0) return;

    if (dp1_f != nullptr)
        delete [] dp1_f;

    delete [] dp1;
    delete [] p1_f;
    delete [] p1;
}
void MomentumGrid::DeallocateP2() {
    if (np2 == 0) return;

    if (dp2_f != nullptr)
        delete [] dp2_f;

    delete [] dp2;
    delete [] p2_f;
    delete [] p2;
}

void MomentumGrid::DeallocatePAndXi0() {
    if (xi0 == nullptr)
        return;
        
    delete [] p;
    delete [] p_f1;
    delete [] p_f2;
	delete [] gamma;
	delete [] gamma_f1;
	delete [] gamma_f2;
    delete [] xi0;
    delete [] xi0_f1;
    delete [] xi0_f2;
 }



/*******************************
 * PUBLIC METHODS              *
 *******************************/
/**
 * Re-build this momentum grid.
 *
 * t:     Time to build this grid for.
 * ri:    Index of point on the radial grid to which this momentum
 *        grid corresponds.
 * rGrid: Radial grid on which this momentum grid lives.
 */
bool MomentumGrid::Rebuild(const real_t t, const len_t ri, const RadialGrid *rGrid) {
    return this->generator->Rebuild(t, ri, this, rGrid);
}


