/**
 * Implementation of the internal boundary condition on the xi = +-1 boundaries.
 * Due to symmetry, the flux across this boundary must vanish, and hence this
 * boundary condition is trivial and doesn't need to set any matrix elements.
 */

#include "FVM/Grid/PXiGrid/XiInternalBoundaryCondition.hpp"

using namespace TQS::FVM::BC;


/**
 * Function called when coefficients for operators should be
 * re-built. Since this boundary condition always the same
 * (zero flux across xi = +-1), we don't need to do anything
 * special here.
 */
bool XiInternalBoundaryCondition::Rebuild(const real_t) { return false; }

/**
 * Set the matrix elements corresponding to this boundary
 * condition. Since the condition is trivial (zero flux across xi = +-1)
 * we don't need to set any matrix element explicitly.
 */
void XiInternalBoundaryCondition::SetMatrixElements(Matrix*) { }

