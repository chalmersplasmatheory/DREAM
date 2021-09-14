/**
 * This term implements the pitch-angle diffusion resulting when electrons
 * interact with the magnetic ripple in a tokamak. The diffusion coefficient
 * was first derived in [Kurzan et al, PRL 75 4626 (1995)].
 */

#include <algorithm>
#include "DREAM/Constants.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/Kinetic/RipplePitchScattering.hpp"
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
 * nCoils:     Number of toroidal field coils.
 * deltaCoils: Distance between toroidal field coils.
 * nModes:     Number of perturbation modes given.
 * m:          Poloidal mode
 * n:          Toroidal mode numbers.
 * dB_B:       Magnetic perturbation levels.
 */
RipplePitchScattering::RipplePitchScattering(
    FVM::Grid *grid, enum OptionConstants::eqterm_ripple_mode mode,
    enum OptionConstants::momentumgrid_type mgtype,
    const len_t nCoils, const len_t nModes, const int_t *m,
    const int_t *n, DREAM::MultiInterpolator1D *dB_B
) : RipplePitchScattering(
    grid, mode, mgtype, (2*M_PI*grid->GetRadialGrid()->GetR0() / nCoils), nModes, m, n, dB_B
) {
    SetName("RipplePitchScattering");

    if (grid->GetRadialGrid()->GetR0() == std::numeric_limits<real_t>::infinity())
        throw DREAMException(
            "PitchScattermTerm: Please use the parameter 'deltaCoils' instead of "
            "'nCoils' in cylindrical mode."
        );
}

RipplePitchScattering::RipplePitchScattering(
    FVM::Grid *grid, enum OptionConstants::eqterm_ripple_mode mode,
    enum OptionConstants::momentumgrid_type mgtype,
    const real_t deltaCoils, const len_t nModes, const int_t *m,
    const int_t *n, DREAM::MultiInterpolator1D *dB_B
) : DiffusionTerm(grid), mode(mode), deltaCoils(deltaCoils), nModes(nModes), m(m), n(n), dB_B(dB_B) {
    
    SetName("RipplePitchScattering");

    if (mgtype != OptionConstants::MOMENTUMGRID_TYPE_PXI)
        throw DREAMException("PitchScatterTerm: Only p-xi grids are supported.");

    Allocate();
}

/**
 * Destructor.
 */
RipplePitchScattering::~RipplePitchScattering() {
    delete [] this->m;
    delete [] this->n;
    delete [] this->p_mn[0];
    delete [] this->p_mn;

    delete this->dB_B;
}

/**
 * Allocate memory used by this term.
 */
void RipplePitchScattering::Allocate() {
    const len_t nr = this->grid->GetNr();

    this->p_mn = new real_t*[this->nModes];
    this->p_mn[0] = new real_t[this->nModes * nr];
    
    for (len_t k = 1; k < this->nModes; k++)
        this->p_mn[k] = this->p_mn[k-1] + nr;

    CalculateResonantMomentum();
}

/**
 * Evaluate the resonant momentum 'p_mn' for all perturbation modes given
 * to this equation term.
 */
void RipplePitchScattering::CalculateResonantMomentum() {
    const real_t
        c  = Constants::c,
        e  = Constants::ec,
        me = Constants::me;

    for (len_t k = 0; k < this->nModes; k++) {
        for (len_t ir = 0; ir < nr; ir++) {
            const real_t G_R0 = this->grid->GetRadialGrid()->GetBTorG(ir);

            if (INCLUDE_POLOIDAL_FIELD) {
                throw NotImplementedException(
                    "RipplePitchScattering: Cannot take the poloidal field component into account for resonant momentum yet."
                );
            } else {
                p_mn[k][ir] = e*G_R0*deltaCoils / (me*c*n[k]*2*M_PI);
            }
        }
    }
}

/**
 * Called when the grid is rebuilt.
 */
bool RipplePitchScattering::GridRebuilt() {
    delete [] this->p_mn[0];
    delete [] this->p_mn;

    Allocate();

    return true;
}

/**
 * Rebuild the diffusion coefficient.
 *
 * t: Simulation time.
 */
void RipplePitchScattering::Rebuild(const real_t t, const real_t, FVM::UnknownQuantityHandler*) {
    const len_t nr = this->grid->GetNr();

    const real_t e  = Constants::ec;
    const real_t me = Constants::me;

    for (len_t k = 0; k < this->nModes; k++) {
        const real_t *dB_mn_B =this->dB_B->Eval(k, t);

        for (len_t ir = 0; ir < nr; ir++) {
            auto mg = this->grid->GetMomentumGrid(ir);
            const len_t np    = mg->GetNp1();
            const len_t nxi   = mg->GetNp2();
            const real_t *p   = mg->GetP1();
            const real_t *p_f = mg->GetP1_f();
            const real_t *dp  = mg->GetDp1();
            const real_t *xi0 = mg->GetP2_f();

            // Magnetic field strength
            //const real_t B = this->grid->GetRadialGrid()->GetFSA_B(ir);
            const real_t B = this->grid->GetRadialGrid()->GetBmin(ir);

            for (len_t j_f = 0; j_f < nxi+1; j_f++){
                const real_t absxi = min(fabs(xi0[j_f]),1.0); // handles precision errors in fabs
                for (len_t i = 0; i < np; i++) {
                    const real_t ppar  = p[i]*absxi;
                    const real_t p2    = p[i]*p[i];
                    const real_t gamma = sqrt(1+p2);
                    
                    real_t Hmn = 0.0;
                    if(mode == OptionConstants::EQTERM_RIPPLE_MODE_BOX){
                        // Solve for the p width of the resonant region
                        const real_t a = sqrt( dB_mn_B[ir] * absxi * sqrt(1-absxi*absxi) );
                        const real_t p_resonance_lo = p_mn[k][ir] / (absxi+a); // lower limit
                        const real_t p_resonance_hi = p_mn[k][ir] / (absxi-a); // upper limit
                        real_t dpBar = min(p_f[i+1], p_resonance_hi ) - max(p_f[i], p_resonance_lo);
                        // Is resonant momentum interval outside of current cell?
                        if (dpBar <= 0)
                            continue;
                        Hmn     = dpBar / (dp[i]*(p_resonance_hi-p_resonance_lo));
                    } else if (mode == OptionConstants::EQTERM_RIPPLE_MODE_GAUSSIAN){
                        real_t p0 = p_mn[k][ir] / absxi;
                        real_t deltaP = p0 * sqrt( dB_mn_B[ir]*sqrt(1-absxi*absxi) );
                        if(deltaP)
                            Hmn = 1.0 / (absxi*dp[i]) * ( gsl_sf_erf( (p_f[i+1]-p0)/deltaP ) - gsl_sf_erf( (p_f[i]-p0)/deltaP ) );
                        else if ( p_f[i+1]>p0 && p_f[i]<=p0 )
                            Hmn = 2.0 / (absxi*dp[i]);
                    }
                    real_t betapar = ppar/gamma;
                    real_t Dperp   = M_PI/32.0 * e*B/me * betapar * dB_mn_B[ir]*dB_mn_B[ir] * Hmn;

                    D22(ir, i, j_f) += (1-xi0[j_f]*xi0[j_f])/p2 * Dperp;
                }
            }
        }
    }
}
