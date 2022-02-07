#include "DREAM/Equations/Kinetic/BraamsKarneyLOperator.hpp"

using namespace DREAM;

static bool DoIntegralBC(len_t i, len_t j, len_t np1, len_t np2) {
    (void)j, (void)np2;
    return i == np1 - 1;
}

BraamsKarneyLOperator::BraamsKarneyLOperator(FVM::Grid *grid, real_t a)
    : EquationTerm(grid), a(a) {
    GridRebuilt();
}

len_t BraamsKarneyLOperator::GetNumberOfNonZerosPerRow() const {
    return 6;
}

bool BraamsKarneyLOperator::SetJacobianBlock(const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t *) {
    if (derivId == unknId) {
        this->SetMatrixElements(jac, nullptr);
        return true;
    }

    return false;
}

template<typename F1>
void BraamsKarneyLOperator::SetElements(F1 X) {
    const len_t nr  = grid->GetNr();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const real_t
            *dp = mg->GetDp1(),
            *dxi = mg->GetDp2(),
            *dp_f = mg->GetDp1_f(),
            *dxi_f = mg->GetDp2_f();
        const real_t
            *p = mg->GetP1(),
            *xi = mg->GetP2(),
            *gamma = mg->GetGamma();

        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            for (len_t j = 0; j < np2; j++) {
                len_t idx = offset + j * np1 + i;

                if (DoIntegralBC(i, j, np1, np2)) {
                    X(idx, idx, -1);

                    continue;
                }

                auto setCoeffAbs = [&] (int i, int j, real_t val) {
                                       len_t offset_idx = offset + j * np1 + i;
                                       X(idx, offset_idx, val);
                                   };

                auto setCoeff = [&] (int i_offset, int j_offset, real_t val) {
                                    setCoeffAbs(i + i_offset, j + j_offset, val);
                                };

                const real_t dp2C = gamma[i] * gamma[i] / dp[i];
                const real_t dxi2C = (1 - xi[j] * xi[j]) / (p[i]*p[i] * dxi[j]);
                const real_t dpC = (2 / p[i] + 3*p[i]);
                const real_t dxiC =-2*xi[j]/(p[i]*p[i]);

                if (i == 0 || j == 0 || j == np2 - 1) {
                    if (i == 0) {
                        setCoeff(1, 0, dp2C / dp_f[0]);
                        setCoeff(0, 0, -dp2C / dp_f[0]);

                        setCoeff(1, 0, dpC / dp_f[0]);
                        setCoeff(0, 0, -dpC / dp_f[0]);
                    } else if (j == 0) {
                        setCoeff(0, +1, dxi2C / dxi_f[0]);
                        setCoeff(0, 0, -dxi2C / dxi_f[0]);

                        setCoeff(0, +1, dxiC / (2 * dxi_f[0]));
                        setCoeff(0, 0, -dxiC / (2 * dxi_f[0]));
                    } else if (j == np2 - 1) {
                        setCoeff(0, -1, dxi2C / dxi_f[j-1]);
                        setCoeff(0, 0, -dxi2C / dxi_f[j-1]);

                        setCoeff(0, 0, dxiC / (2 * dxi_f[j - 1]));
                        setCoeff(0, -1, -dxiC / (2 * dxi_f[j - 1]));
                    }
                } else {
                    // d^2 psi / d p^2
                    setCoeff(+1, 0, dp2C / dp_f[i]);
                    setCoeff(-1, 0, dp2C / dp_f[i-1]);
                    setCoeff( 0, 0, -dp2C * (dp_f[i]+dp_f[i-1]) / (dp_f[i]*dp_f[i-1]));

                    // d^2 psi / d xi^2
                    setCoeff(0, +1, dxi2C / dxi_f[j]);
                    setCoeff(0, -1, dxi2C / dxi_f[j-1]);
                    setCoeff(0,  0, -dxi2C * (dxi_f[j]+dxi_f[j-1]) / (dxi_f[j]*dxi_f[j-1]));

                    // d psi / d p
                    setCoeff(+1, 0, dpC / (dp_f[i] + dp_f[i-1]));
                    setCoeff(-1, 0,-dpC / (dp_f[i] + dp_f[i-1]));

                    // d psi / d xi
                    setCoeff(0, +1, dxiC / (dxi_f[j] + dxi_f[j-1]));
                    setCoeff(0, -1,-dxiC / (dxi_f[j] + dxi_f[j-1]));

                    // psi
                    setCoeff(0, 0, 1 - a*a);
                }
            }
        }

        offset += np1 * np2;
    }
}

void BraamsKarneyLOperator::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    SetElements([&](auto idx, auto idx_p, auto V) {
                    mat->SetElement(idx, idx_p, V);
                });
}

void BraamsKarneyLOperator::SetVectorElements(real_t *vec, const real_t *f) {
    SetElements([&](auto idx, auto idx_p, auto V) {
                    vec[idx] += f[idx_p] * V;
                });
}
