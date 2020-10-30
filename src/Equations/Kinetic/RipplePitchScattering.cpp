/**
 * This term implements the pitch-angle diffusion resulting when electrons
 * interact with the magnetic ripple in a tokamak. The diffusion coefficient
 * was first derived in [Kurzan et al, PRL 75 4626 (1995)].
 */

#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Kinetic/RipplePitchScattering.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Interpolator1D.hpp"


using namespace DREAM;


/**
 * Constructor.
 *
 * grid:   Kinetic grid to which this term is applied.
 * nCoils: Number of toroidal field coils.
 * nModes: Number of perturbation modes given.
 * m:      Poloidal mode
 * n:      Toroidal mode numbers.
 * dB_B:   Magnetic perturbation levels.
 */
RipplePitchScattering::RipplePitchScattering(
    FVM::Grid *grid, const len_t nCoils, const len_t nModes,
    const len_t *m, const len_t *n, FVM::Interpolator1D **dB_B
) : DiffusionTerm(grid), nCoils(nCoils), nModes(nModes), m(m), n(n), dB_B(dB_B) {

    Allocate();
}

/**
 * Destructor.
 */
RipplePitchScattering::~RipplePitchScattering() {
    delete [] this->m;
    delete [] this->n;
    delete [] this->p_mn;

    for (len_t i = 0; i < this->nModes; i++)
        delete this->dB_B[i];

    delete [] this->dB_B;
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
            const real_t G = this->grid->GetRadialGrid()->GetBTorG(ir);

            if (INCLUDE_POLOIDAL_FIELD) {
                throw NotImplementedException(
                    "RipplePitchScattering: Cannot take the poloidal field component into account for resonant momentum yet."
                );
            } else {
                p_mn[k][ir] = e*G / (me*c*n[k]*nCoils);
            }
        }
    }
}

/**
 * Called when the grid is rebuilt.
 */
bool RipplePitchScattering::GridRebuilt() {
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

    const real_t c  = Constants::c;
    const real_t e  = Constants::ec;
    const real_t me = Constants::me;

    for (len_t k = 0; k < this->nModes; k++) {
        const real_t *dB_mn_B =this->dB_B[k]->Eval(t);

        for (len_t ir = 0; ir < nr; ir++) {
            auto mg = this->grid->GetMomentumGrid(ir);
            const len_t np1   = mg->GetNp1();
            const len_t np2   = mg->GetNp2();
            const real_t *p   = mg->GetP1();
            const real_t *xi0 = mg->GetP2_f();

            // Magnetic field strength
            //const real_t B = this->grid->GetRadialGrid()->GetFSA_B(ir);
            const real_t B = this->grid->GetRadialGrid()->GetBmin(ir);

            for (len_t j_f = 0; j_f < np2+1; j_f++) {
                for (len_t i = 0; i < np1; i++) {
                    const real_t ppar  = p[i]*xi0[j_f];
                    const real_t p2    = p[i]*p[i];
                    const real_t gamma = sqrt(1+p2);

                    // Resonant momentum
                    const real_t delta_p_mn = sqrt(dB_mn_B[ir] * ppar * sqrt(p2-ppar*ppar));

                    if (fabs(ppar - p_mn[k][ir]) >= delta_p_mn)
                        continue;

                    real_t Hmn   = 1/delta_p_mn;
                    real_t vpar  = c*ppar/gamma;
                    real_t Dperp = M_PI/32.0 * e*B/me * vpar * dB_mn_B[ir]*dB_mn_B[ir] * Hmn;

                    D22(ir, i, j_f) += (1-xi0[j_f]*xi0[j_f])/p2 * Dperp;
                }
            }
        }
    }
}

