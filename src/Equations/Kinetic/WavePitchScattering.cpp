/**
 * This term implements the pitch-angle diffusion resulting when electrons
 * interact with the magnetic ripple in a tokamak. The diffusion coefficient
 * was first derived in [Kurzan et al, PRL 75 4626 (1995)].
 */

#include <algorithm>
#include "DREAM/Constants.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/WavePitchScattering.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Interpolator1D.hpp"
#include <gsl/gsl_sf_erf.h>
using namespace DREAM;
using namespace std;


/**
 * Constructor.
 *
 * grid:       Kinetic grid to which this term is applied.
 * ppar_res:        Resonant momentum (vs minor radius)
 * Delta_ppar_res:       Resonance width (vs minor radius)
 * Dxx_int:     integrated pitch angle diffusion term (vs minor radius)
 * t_start:     starting time for wave scattering
 * t_end:       end time for wave scattering
 */
WavePitchScattering::WavePitchScattering(
    FVM::Grid *grid, enum OptionConstants::eqterm_wave_mode mode,
    enum OptionConstants::momentumgrid_type mgtype,
    DREAM::MultiInterpolator1D *ppar_res, DREAM::MultiInterpolator1D *Delta_ppar_res, DREAM::MultiInterpolator1D *Dxx_int
) : WavePitchScattering(
    grid, mode, mgtype, ppar_res, Delta_ppar_res, Dxx_int
) {
    SetName("WavePitchScattering");
}

WavePitchScattering::WavePitchScattering(
    FVM::Grid *grid, enum OptionConstants::eqterm_wave_mode mode,
    enum OptionConstants::momentumgrid_type mgtype,
    DREAM::MultiInterpolator1D *ppar_res, DREAM::MultiInterpolator1D *Delta_ppar_res, DREAM::MultiInterpolator1D *Dxx_int
) : DiffusionTerm(grid), mode(mode), ppar_res(ppar_res), Delta_ppar_res(Delta_ppar_res), Dxx_int(Dxx_int)
{
    
    SetName("WavePitchScattering");

    if (mgtype != OptionConstants::MOMENTUMGRID_TYPE_PXI){
            throw DREAMException("PitchScatterTerm: Only p-xi grids are supported.");
        }
        
    Allocate();
}

/**
 * Allocate memory used by this term.
 */
void WavePitchScattering::Allocate() {
    // nothing to be allocated
}

/**
 * Destructor.
 */
WavePitchScattering::~WavePitchScattering() {
    delete this->ppar_res;
    delete this->Delta_ppar_res;
    delete this->Dxx_int;
}

/**
 * Called when the grid is rebuilt.
 */
bool WavePitchScattering::GridRebuilt() {
    // nothing to be done
    Allocate();
    return true;
}

/**
 * Rebuild the diffusion coefficient at specified simulation time.
 *
 * t: Simulation time.
 */
void WavePitchScattering::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*) {
    const len_t nr = this->grid->GetNr();
    
    // evaluate wave quantities at this time
    const real_t *ppar_res_t =this->ppar_res->Eval(1, t);
    const real_t *Delta_ppar_res_t =this->Delta_ppar_res->Eval(1, t);
    const real_t *Dxx_int_t =this->Dxx_int->Eval(1, t);
    
    // loop over radial coordinate
    for (len_t ir = 0; ir < nr; ir++) {
        auto mg = this->grid->GetMomentumGrid(ir);
        const len_t np    = mg->GetNp1();
        const len_t nxi   = mg->GetNp2();
        const real_t *p   = mg->GetP1();
        const real_t *xi0 = mg->GetP2_f();
            
        // loop over pitch angle grid
        for (len_t j_f = 0; j_f < nxi+1; j_f++){
            const real_t absxi = min(fabs(xi0[j_f]),1.0); // handles precision errors in fabs
            
            // loop over momentum grid
            for (len_t i = 0; i < np; i++) {

                if (mode == OptionConstants::EQTERM_WAVE_MODE_GAUSSIAN){
                    if (absxi != 0) {
                        // calculate resonant momentum and total diffusion strength at this pitch angle
                        real_t p_res0 = ppar_res_t[ir] / absxi;
                        real_t Delta_p_res0 = Delta_ppar_res_t[ir] / absxi;
                        real_t Dxx_int0 = Dxx_int_t[ir];
                        
                        // add to diffusion term
                        const real_t exp_nom = (p[i]-p_res0)*(p[i]-p_res0);
                        const real_t exp_denom = 2*Delta_p_res0*Delta_p_res0;
                        D22(ir, i, j_f) += Dxx_int0 * exp(-exp_nom/exp_denom);
                    }
                }
            }
        }
    }
}
