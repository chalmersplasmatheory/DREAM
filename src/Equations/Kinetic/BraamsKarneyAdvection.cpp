#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyAdvection.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

void BraamsKarneyAdvection::DCoeffs::init(len_t nr, len_t *n1, len_t *n2) {
    this->df1 = new real_t*[nr];
    this->df2 = new real_t*[nr];

    len_t
        nElements_f1 = 0,
        nElements_f2 = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
    }

	this->df1[0] = new real_t[nElements_f1];
	this->df2[0] = new real_t[nElements_f2];

	for (len_t i = 1; i < nr; i++) {
		this->df1[i] = this->df1[i-1] + ((n1[i-1]+1)*n2[i-1]);
		this->df2[i] = this->df2[i-1] + (n1[i-1]*(n2[i-1]+1));
	}
}

BraamsKarneyAdvection::DCoeffs::~DCoeffs() {
	if (df1 != nullptr) {
        delete [] df1[0];
        delete [] df1;
    }
    if (df2 != nullptr) {
        delete [] df2[0];
        delete [] df2;
    }
}

/**
 * Constructor.
 */
BraamsKarneyAdvection::BraamsKarneyAdvection(FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, CoulombLogarithm *lnLambda, len_t id_pi_0, len_t id_pi_1)
    : FVM::AdvectionTerm(grid, unknowns), lnLambda(lnLambda), id_pi_0(id_pi_0), id_pi_1(id_pi_1) {

    SetName("BraamsKarneyAdvection");

    AddUnknownForJacobian(unknowns, id_pi_0);
    AddUnknownForJacobian(unknowns, id_pi_1);

	for (int i = 0; i < dfs_n1; i++)
		for (int j = 0; j < dfs_n2; j++)
		dfs[j][i].init(nr, n1, n2);
}

static real_t Gamma(const real_t p) {
    return sqrt(1 + p * p);
}

const real_t alphaBarOverLnLambda = 3.75927427447396469133e-19; // ec^4 / (eps0^2 * me^2 * c^3)

template<typename T1>
void BraamsKarneyAdvection::SetCoefficients(T1 psi, real_t **f1, real_t **f2) {
	len_t r_offset = 0;

	for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);

        const len_t np1 = n1[ir];
        const len_t np2 = n2[ir];

		const real_t
			*xi_f = mg->GetP2_f(),
			*p = mg->GetP1(),
			*p_f = mg->GetP1_f(),
			*gamma = mg->GetGamma(),
			*dp_f = mg->GetDp1(),
			*dxi_f = mg->GetDp2()
            ;

        //real_t alphabar = alphaBarOverLnLambda * lnLambda->evaluateLnLambdaC(ir);

        CoulombLogarithm::collqty_settings settings;
        settings.lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
        auto alphabar = [&] (real_t p) {
                            return alphaBarOverLnLambda * lnLambda->evaluateAtP(ir, p, &settings);
                        };


        auto callPsi = [&](len_t i, len_t j, int di, int dj) {
                           len_t idx = r_offset + (i + di) + np1 * (j + dj);
                           return psi(idx, di , dj);
                       };

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1+1; i++) {
				if (i == 0) {
					// Zero.
				} else if (i == np1) {
					// Same as i = np1 - 1.
					real_t pi_combo_low = callPsi(i, j, -2, 0);
					real_t pi_combo_high = callPsi(i, j, -1, 0);
					F1(ir, i, j, f1) = -alphabar(p_f[i]) * Gamma(p_f[i - 1]) * (pi_combo_high - pi_combo_low) / dp_f[i - 1];
				} else { // Not on a boundary.
					real_t pi_combo_low = callPsi(i, j, -1, 0);
					real_t pi_combo_high = callPsi(i, j, 0, 0);
					F1(ir, i, j, f1) = -alphabar(p_f[i]) * Gamma(p_f[i]) * (pi_combo_high - pi_combo_low) / dp_f[i];
				}
			}
		}

		for (len_t j = 0; j < np2+1; j++)
			for (len_t i = 0; i < np1; i++) {
				if (j == 0 || j == np2) {
					// Zero.
				} else {
					real_t pi_combo_low = callPsi(i, j, 0, -1);
					real_t pi_combo_high = callPsi(i, j, 0, 0);

                    F2(ir, i, j, f2) = alphabar(p[i]) * sqrt(1 - xi_f[j] * xi_f[j]) / (p[i] * gamma[i])
                        * (pi_combo_high - pi_combo_low) / dxi_f[j];
                }
			}

		r_offset += np1 * np2;
	}
}

/**
 * Build the coefficients of this advection (or diffusion) term.
 */
void BraamsKarneyAdvection::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *x){
    real_t *pi_0 = x->GetUnknownData(id_pi_0);
    real_t *pi_1 = x->GetUnknownData(id_pi_1);

    SetCoefficients(
        [&] (len_t idx, int, int) {
            return pi_0[idx] - 2 * pi_1[idx];
        }, this->f1, this->f2);
}


// Set jacobian of the advection coefficients for this advection term
void BraamsKarneyAdvection::SetPartialAdvectionTerm(len_t derivId, len_t /*nMultiples*/){
    ResetDifferentiationCoefficients();

	for (int i = 0; i < dfs_n1; i++) {
		for (int j = 0; j < dfs_n2; j++) {
			real_t psi_coeff = 1;
			if (derivId == id_pi_1) {
				psi_coeff = -2;
			} else if (derivId != id_pi_0) {
				std::abort();
			}
			SetCoefficients(
				[&] (len_t, int di, int dj) {
					return (di == n1_offsets[i] && dj == n2_offsets[j]) ? psi_coeff : 0;
				}, this->dfs[j][i].df1, this->dfs[j][i].df2);
		}
	}
}

// Overridden
void BraamsKarneyAdvection::SetPartialJacobianContribution(int_t diagonalOffset, jacobian_interp_mode, len_t, FVM::Matrix *jac, const real_t *x, bool){
	if (diagonalOffset)
		return;

	for (int k1 = 0; k1 < dfs_n1; k1++) {
		for (int k2 = 0; k2 < dfs_n2; k2++) {
			int n1_offset = n1_offsets[k1];
			int n2_offset = n2_offsets[k2];

			// if (n2_offset) // Only set main diagonal.
			// 	continue;

			ResetJacobianColumn();

			auto setVectorElements = [&] (
				real_t *vec, const real_t *x,
				int k1, int k2
				) {
#define f(K,I,J,V) vec[offset+j*np1+i] += (V)*x[offset+((K)-ir)*np2*np1 + (J)*np1 + (I)]
#   include "BraamsKarneyAdvection.set.cpp"
#undef f
			};

			setVectorElements(JacobianColumn, x, k1, k2);

			len_t offset = 0;
			for(len_t ir=0; ir<nr; ir++){
				for (len_t j = 0; j < n2[ir]; j++)
					for (len_t i = 0; i < n1[ir]; i++) {
						if ((int)i + n1_offset < 0 || (int)i + n1_offset >= (int)n1[ir] ||
							(int)j + n2_offset < 0 || (int)j + n2_offset >= (int)n2[ir])
							continue;
						jac->SetElement(offset + n1[ir]*j + i, offset + n1[ir]*(j + n2_offset) + i + n1_offset, JacobianColumn[offset + n1[ir]*j + i]);
					}
				offset += n1[ir]*n2[ir];
			}
		}
	}
}
