/**
 * Routines for evaluating the equilibrium charge state distribution
 * of ions.
 */

#include <vector>
#include "DREAM/ADAS.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;
using namespace std;


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
		return (std::abs(x1) < std::numeric_limits<real_t>::min() / std::abs(x2));
	} else // x2 == 1
		return false;
}
#endif

/**
 * Initialize the ion densities for the specified species according to
 * coronal equilibrium.
 *
 * ih:            Ion handler.
 * adas:          ADAS rate handler.
 * equil_indices: Ion indices for the ions to be initialized in equilibrium.
 * nfree0:        Vector containing background free electron density as function
 *                of radius.
 * Ni:            Total ion densities.
 * T:             Electron temperature.
 * Nr:            Number of radial grid points.
 *
 * RETURNS
 * ni:            List of ion charge state and radial density distributions.
 */
void SimulationGenerator::EvaluateIonEquilibrium(
	IonHandler *ih, ADAS *adas, vector<len_t>& equil_indices,
	const real_t *nfree0, const real_t *Ni, const real_t *T,
	len_t Nr, vector<real_t*>& ni
) {
	len_t Nequil = equil_indices.size();
	vector<len_t> Z(Nequil);
	
	// Allocate ADAS rates
	vector<real_t*> I(Nequil);
	vector<real_t*> R(Nequil);
	for (len_t i = 0; i < Nequil; i++) {
		Z[i] = ih->GetZ(equil_indices[i]);
		I[i] = new real_t[Z[i]];
		R[i] = new real_t[Z[i]];
	}

	// Initial guess for nfree is that all atoms are fully ionized
	real_t eps = std::numeric_limits<real_t>::epsilon();
	const len_t MAX_ITER=50;
	len_t nmaxiter = 0;
	for (len_t ir = 0; ir < Nr; ir++) {
		real_t nfree = nfree0[ir], dnfree;

		// Assume equilibrium ions are singly ionized
		for (len_t & iZ : equil_indices)
			nfree += Ni[iZ*Nr+ir];

		len_t niter = 0;
		do {
			for (len_t iZ = 0; iZ < Nequil; iZ++) {
				const real_t *Ntot = Ni+equil_indices[iZ]*Nr;
				real_t *II = I[iZ];
				real_t *RR = R[iZ];

				EvaluateADASRates(adas, Z[iZ], nfree, T[ir], II, RR);

				for (len_t l = 0; l < Z[iZ]+1; l++) {
					real_t s = 0;
					if (l > 0) {
						for (len_t j = 0; j < l; j++) {
							real_t p = 1;
							for (len_t k = j+1; k <= l; k++) {
#if !defined(NDEBUG) && defined(__linux__)
								if (OUFLOW(p, RR[k]/II[k-1]))
									break;
#endif
								p *= RR[k] / II[k-1];
							}

							s += p;
						}
					}

					if (l < Z[iZ]) {
						for (len_t j = l+1; j <= Z[iZ]; j++) {
							real_t p = 1;
							for (len_t k = l; k < j; k++) {
#if !defined(NDEBUG) && defined(__linux__)
								if (OUFLOW(p, II[k]/RR[k+1]))
									break;
#endif
								p *= II[k] / RR[k+1];
							}

							s += p;
						}
					}

					ni[iZ][l*Nr+ir] = Ntot[ir] / (1 + s);
				}
			}
			
			// Evaluate free electron density
			real_t nfree1 = nfree0[ir];
			for (len_t iZ = 0; iZ < Nequil; iZ++)
				for (len_t j = 1; j < Z[iZ]+1; j++)
					nfree1 += j * ni[iZ][j*Nr+ir];

			dnfree = std::abs(nfree1-nfree);
			nfree = nfree1;

			niter++;
		} while (dnfree > nfree*eps && niter < MAX_ITER);

		if (niter > nmaxiter)
			nmaxiter = niter;
	}

	//printf("Number of iterations: " LEN_T_PRINTF_FMT "\n", nmaxiter);
}

