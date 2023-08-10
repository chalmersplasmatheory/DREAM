/**
 * Routines for evaluating the equilibrium charge state distribution
 * of ions.
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;


/**
 * Evaluate the ADAS rates at the given density and temperature.
 */
void SimulationGenerator::EvaluateADASRates(
	ADAS *adas, const len_t Z, const real_t n, const real_t T,
	real_t *I, real_t *R
) {
	ADASRateInterpolator *acd = adas->GetACD(Z);
	ADASRateInterpolator *scd = adas->GetSCD(Z);

	for (len_t j = 0; j < Z+1; j++) {
		I[j] = scd->Eval(j, n, T);
		R[j] = acd->Eval(j, n, T);
	}
}

#if !defined(NDEBUG) && defined(__linux__)
/**
 * Helper function to check if the product of two floating-point values will
 * over/underflow. This is required when compiling DREAM in debug mode, since
 * floating-point exceptions are then enabled, preventing the desired "inf"
 * to be produced at overflow (or 0 at underflow).
 */
bool OUFLOW(real_t x1, real_t x2) {
	if (x2 > 1) {
		return (std::abs(x1) > std::numeric_limits<real_t>::max() / std::abs(x2));
	} else if (x2 < 1) {
		return (std::abs(x1) > std::numeric_limits<real_t>::min() / std::abs(x2));
	} else // x2 == 1
		return false;
}
#endif

/**
 * Initialize the ion densities for the specified species according to
 * coronal equilibrium.
 *
 * ih:   Ion handler.
 * adas: ADAS rate handler.
 * idx:  Index of ion species in IonHandler.
 * Ni:   Total ion density.
 * T:    Electron temperature.
 * Nr:   Number of radial grid points.
 *
 * RETURNS
 * ni:   Ion charge state and radial density distribution.
 */
void SimulationGenerator::EvaluateIonEquilibrium(
	IonHandler *ih, ADAS *adas, len_t idx,
	const real_t *Ni, const real_t *T, len_t Nr, real_t *ni
) {
	len_t Z = ih->GetZ(idx);
	
	// Evaluate ADAS rates
	real_t *I = new real_t[Z];
	real_t *R = new real_t[Z];
	
	// Initial guess for nfree is that all atoms are fully ionized
	real_t eps = std::numeric_limits<real_t>::epsilon();
	const len_t MAX_ITER=20;
	len_t nmaxiter = 0;
	for (len_t ir = 0; ir < Nr; ir++) {
		real_t nfree = Ni[ir], dnfree;
		len_t niter = 0;
		do {
			EvaluateADASRates(adas, Z, nfree, T[ir], I, R);

			for (len_t l = 0; l < Z+1; l++) {
				real_t s = 0;
				if (l > 0) {
					for (len_t j = 0; j < l; j++) {
						real_t p = 1;
						for (len_t k = j+1; k <= l; k++) {
#if !defined(NDEBUG) && defined(__linux__)
							if (OUFLOW(p, R[k]/I[k-1]))
								break;
#endif
							p *= R[k] / I[k-1];
						}

						s += p;
					}
				}

				if (l < Z) {
					for (len_t j = l+1; j <= Z; j++) {
						real_t p = 1;
						for (len_t k = l; k < j; k++) {
#if !defined(NDEBUG) && defined(__linux__)
							if (OUFLOW(p, I[k]/R[k+1]))
								break;
#endif
							p *= I[k] / R[k+1];
						}

						s += p;
					}
				}

				ni[l*Nr+ir] = Ni[ir] / (1 + s);
			}
			
			// Evaluate free electron density
			real_t nfree1 = 0;
			for (len_t j = 1; j < Z+1; j++)
				nfree1 += j * ni[j*Nr+ir];

			dnfree = std::abs(nfree1-nfree);
			nfree = nfree1;

			niter++;
		} while (dnfree > nfree*eps && niter < MAX_ITER);

		if (niter > nmaxiter)
			nmaxiter = niter;
	}

	printf("Number of iterations: " LEN_T_PRINTF_FMT "\n", nmaxiter);
}

