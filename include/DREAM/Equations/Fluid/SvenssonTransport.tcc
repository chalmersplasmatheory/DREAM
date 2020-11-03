/**
 * Implementation of a diffusive transport term which can be applied to both
 * kinetic and fluid grids, and which allows one to prescribe the diffusion
 * coefficient in time and phase space.
 */

#include <type_traits>
#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"

/**
 * Constructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::SvenssonTransport(
    DREAM::FVM::Grid *grid,
    const len_t nr, const len_t np, const real_t pStar,
    const real_t **coeffA, const real_t **coeffD,
    const real_t *r, const real_t *p,
    DREAM::FVM::UnknownQuantityHandler *unknowns,
    DREAM::RunawayFluid *REFluid,
    bool allocCoefficients
) : T(grid, allocCoefficients),
    nr(nr), np(np), pStar(pStar),
    coeffA(coeffA), coeffD(coeffD), r(r), p(p),
    unknowns(unknowns), REFluid(REFluid)
{
//protected:
    this->EID = this->unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD); 
}

/**
 * Destructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::~SvenssonTransport() {
    // YYY I guess that something should be in here?
    // if (this->coeffD != nullptr)
    //     delete this->coeffD;
    // if (this->coeffA != nullptr)
    //     delete this->coeffA;
}



/**
 * Rebuild this term by evaluating and setting the diffusion
 * coefficient for the next time step.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::Rebuild(
    const real_t, const real_t, DREAM::FVM::UnknownQuantityHandler*
    ) {
    //const real_t *c = this->prescribedCoeff->Eval(t);
    
    const len_t nr = this->grid->GetNr();
    
    
    // Iterate over the radial flux grid...
    for (len_t ir = 0, offset = 0; ir < nr+1; ir++) {
        // Need interpolation from cell grid to flux grid:
        // pBar_f[0]=pBar[0]
        // pBar_f[ir]  = (pBar[ir-1] + pBar[ir] )*.5
        // pBar_f[nr]= extrapolate
        
        // The varaible to be added to 
        real_t pIntCoeff = 0;
        const real_t *integrandArray = this->EvaluateIntegrand(ir);
            for (len_t i = 0; i < this->np; i++) {
                // The actual integration in p
                //pIntCoeff += integrandArray[i+offset];
                pIntCoeff += integrandArray[i]; // I think that the offset 
                                                // is now baked into
                                                // EvaluateIntegrand(ir)
            }
        // Sets the 
        this->_setcoeff(ir, pIntCoeff);
        offset += this->np;
    }
}

