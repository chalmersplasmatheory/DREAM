#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyAdvection.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
BraamsKarneyAdvection::BraamsKarneyAdvection(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns,
											 len_t id_f, len_t id_pi_0, len_t id_pi_1)
    : FVM::AdvectionTerm(g), id_f(id_f), id_pi_0(id_pi_0), id_pi_1(id_pi_1) {

    SetName("BraamsKarneyAdvection");

    AddUnknownForJacobian(unknowns, id_pi_0);
    AddUnknownForJacobian(unknowns, id_pi_1);
}


/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void BraamsKarneyAdvection::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x){
    real_t *pi_0 = x->GetUnknownData(id_pi_0);
    real_t *pi_1 = x->GetUnknownData(id_pi_1);

	len_t r_offset = 0;

	for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);

        const len_t np1 = n1[ir];
        const len_t np2 = n2[ir];

		const real_t *gamma_f1 = mg->GetGamma_f1(),
			*gamma = mg->GetGamma(),
			*p = mg->GetP1(),
			*xi_f = mg->GetP2_f(),
			*dp_f = mg->GetDp1(),
			*dxi_f = mg->GetDp2();

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
				if (i == 0) {
					// Zero.
				} else if (i == np1) {
					// Same as i = np1 - 1.
					real_t pi0_low = pi_0[r_offset + j*np1 + i - 2];
					real_t pi0_high = pi_0[r_offset + j*np1 + i - 1];
					real_t pi1_low = pi_1[r_offset + j*np1 + i - 2];
					real_t pi1_high = pi_1[r_offset + j*np1 + i - 1];
					real_t pi_combo_low = pi0_low - 2 * pi1_low;
					real_t pi_combo_high = pi0_high - 2 * pi1_high;
					F1(ir, i, j) += -gamma_f1[i] * (pi_combo_high - pi_combo_low) / dp_f[i];
				} else { // Not on a boundary.
					real_t pi0_low = pi_0[r_offset + j*np1 + i - 1];
					real_t pi0_high = pi_0[r_offset + j*np1 + i];
					real_t pi1_low = pi_1[r_offset + j*np1 + i - 1];
					real_t pi1_high = pi_1[r_offset + j*np1 + i];
					real_t pi_combo_low = pi0_low - 2 * pi1_low;
					real_t pi_combo_high = pi0_high - 2 * pi1_high;

					F1(ir, i, j) += -gamma_f1[i] * (pi_combo_high - pi_combo_low) / dp_f[i];
				}
			}
		}

		for (len_t j = 0; j < np2+1; j++)
			for (len_t i = 0; i < np1; i++) {
				if (j == 0 || j == np2) {
					// Zero.
				} else {
					real_t pi0_low = pi_0[r_offset + (j - 1)*np1 + i];
					real_t pi0_high = pi_0[r_offset + j*np1 + i];
					real_t pi1_low = pi_1[r_offset + (j - 1)*np1 + i];
					real_t pi1_high = pi_1[r_offset + j*np1 + i];
					real_t pi_combo_low = pi0_low - 2 * pi1_low;
					real_t pi_combo_high = pi0_high - 2 * pi1_high;

					F2(ir, i, j) += sqrt(1 - xi_f[j] * xi_f[j]) / (p[i] * gamma[i])
						* (pi_combo_high - pi_combo_low) / dxi_f[j];
				}
			}

		r_offset += np1 * np2;
	}
}


// Set jacobian of the advection coefficients for this advection term
void BraamsKarneyAdvection::SetPartialAdvectionTerm(len_t derivId, len_t /*nMultiples*/){
    ResetDifferentiationCoefficients();
	len_t r_offset = 0;

	// There are two cases: derivId == id_pi_0 or id_pi_1.
    real_t potential_factor = 0;
    if (derivId == id_pi_0) {
        potential_factor = -2;
    } else if (derivId == id_pi_0) {
        potential_factor = 1;
    }

	for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);

        const len_t np1 = n1[ir];
        const len_t np2 = n2[ir];

		const real_t *gamma_f1 = mg->GetGamma_f1(),
			*gamma = mg->GetGamma(),
			*p = mg->GetP1(),
			*xi_f = mg->GetP2_f(),
			*dp_f = mg->GetDp1(),
			*dxi_f = mg->GetDp2();

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
				if (i == 0) {
					// Zero.
				} else if (i == np1) {
					// Same as i = np1 - 1.
                    dF1(ir, i, j, 0) = -gamma_f1[i] * potential_factor / dp_f[i];
                } else { // Not on a boundary.
					dF1(ir, i, j, 0) = -gamma_f1[i] * potential_factor / dp_f[i];
				}
			}
		}

		for (len_t j = 0; j < np2+1; j++)
			for (len_t i = 0; i < np1; i++) {
				if (j == 0 || j == np2) {
					// Zero.
				} else {
					dF2(ir, i, j, 0) = sqrt(1 - xi_f[j] * xi_f[j]) / (p[i] * gamma[i])
						* potential_factor / dxi_f[j];
				}
			}

		r_offset += np1 * np2;
	}
}



