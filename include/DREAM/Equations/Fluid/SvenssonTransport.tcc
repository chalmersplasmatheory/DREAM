/**
 * Implementation of a diffusive transport term which can be applied to both
 * kinetic and fluid grids, and which allows one to prescribe the diffusion
 * coefficient in time and phase space.
 */

#include <type_traits>
#include "DREAM/Equations/Fluid/SvennsonTransport.hpp"
//#include "FVM/Interpolator1D.hpp"
//#include "FVM/Interpolator3D.hpp"


/**
 * Constructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::SvenssonTransport(
    DREAM::FVM::Grid *grid,
    const len_t nr, const len_t np, const real_t pStar,
    const real_t **coeffA, const real_t **coeffD,
    const real_t *r, const real_t *p,
    DREAM::FVM::UnknownQuantityHandler* unknowns,
    DREAM::RunawayFluid* REFluid,
    bool allocCoefficients
) : T(grid, allocCoefficients),
    nt(nt), nr(nr), np(np),
    coeffA(coeffA), coeffD(coeffD), t(t), r(r), p(p),
    REFluid(REFluid)
{
    len_t this->EID = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD); 
}

/**
 * Destructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::~SvenssonTransport() {
    if (this->prescribedCoeff != nullptr)
        delete this->prescribedCoeff;
    if (this->interpolateddata != nullptr) {
        delete [] this->interpolateddata[0];
        delete [] this->interpolateddata;
    }
}


/**
 * Function called whenever the computational grid has beenx
 * rebuilt. Here, we will need to interpolate the prescribed
 * diffusion coefficient onto the new grid.
 */
// template<typename T>
// bool DREAM::SvenssonTransport<T>::GridRebuilt() {
//     //  Call GridRebuilt() on base class...
//     this->T::GridRebuilt();

//     // Interpolate prescribed coefficient onto new
//     // phase space grid...
//     InterpolateCoefficient();

//     return true;
// }

/**
 * Interpolate the input coefficient onto the phase space
 * grid used by this EquationTerm.
 */
// template<typename T>
// void DREAM::SvenssonTransport<T>::InterpolateCoefficient() {
//     real_t **newdata = new real_t*[nt];
//     // Drr is defined on the radial flux grid so we
//     // XXX assume that all momentum grids are the same
//     // and one momentum grid's worth of cells to the total
//     // number of cells...
// }

/**
 * Rebuild this term by evaluating and setting the diffusion
 * coefficient for the next time step.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::Rebuild(
    const real_t t, const real_t, DREAM::FVM::UnknownQuantityHandler*
    ) {
    //const real_t *c = this->prescribedCoeff->Eval(t);
    
    const len_t nr = this->grid->GetNr();
    
    
    // Iterate over the radial flux grid...
    for (len_t ir = 0, offset = 0; ir < nr+1; ir++) {
	// Need interpolation from cell grid to flux grid:
	// pBar_f[0]=pBar[0]
	// avg [i]  = (pBar[i-1] + pBar[i] )*.5
	// pBar_f[nr]= extrapolate

	// The varaible to be added to 
	real_t pIntCoeff = 0;
	const reat_t* givenCoeff = this->EvaluateIntegrand(ir)
	    for (len_t i = 0; i < this->np; i++) {
		// The actual integration in p
		// pIntCoeff+=... givenCoeff[i+offset]
	    }
	// Sets the 
	this->_setcoeff(ir, pIntCoeff);
	offset += this->np;
    }
}

