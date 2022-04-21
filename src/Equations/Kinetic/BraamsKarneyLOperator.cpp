#include "DREAM/Equations/Kinetic/BraamsKarneyLOperator.hpp"

#include <iostream>

using namespace DREAM;

static bool DoIntegralBC(len_t i, len_t j, len_t np1, len_t np2) {
    (void)j, (void)np2;
	//return false;
    //return true;
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
        // const real_t
        //     *dp = mg->GetDp1(),
        //     *dxi = mg->GetDp2(),
        //     *dxi_f = mg->GetDp2();
        // const real_t
        //     *p = mg->GetP1(),
        //     *xi = mg->GetP2(),
        //     *gamma = mg->GetGamma();


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

        // real_t dxi0 = 2 * (xi[0] + 1);
        // real_t dximax = 2 * (1 - xi[np2 - 1]);

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

                // const real_t dp2C = gamma[i] * gamma[i] / dp[i];
                // const real_t dxi2C = (1 - xi[j] * xi[j]) / (p[i]*p[i] * dxi[j]);
                // const real_t dpC = (2 / p[i] + 3*p[i]);
                // const real_t dxiC =-2*xi[j]/(p[i]*p[i]);
                // const real_t Cpd1 = 2 / p[i] + 3 * p[i];
                // const real_t Cpd2 = gamma[i] * gamma[i];

                // const real_t Cxid1 = -2 * xi[j] / (p[i] * p[i]);
                // const real_t Cxid2 = (1 - xi[j] * xi[j]) / (p[i] * p[i]);
                // // const real_t Cdxi2 = gamma[i] * gamma[i];
                // // const real_t dpC = (2 / p[i] + 3*p[i]);
                // // const real_t dxiC =-2*xi[j]/(p[i]*p[i]);

                // if (j == 0 || j == np2 - 1) {
                //     if (j == 0) {
                //         setCoeff(0, +1, Cxid2 / dxi_f[0]);
                //         setCoeff(0, 0, -Cxid2 / dxi_f[0]);

                //         setCoeff(0, +1, Cxid1 / (2 * dxi_f[0]));
                //         setCoeff(0, 0, -Cxid1 / (2 * dxi_f[0]));
                //     } else if (j == np2 - 1) {
                //         setCoeff(0, 0, Cxid2 / dxi_f[j-1]);
                //         setCoeff(0, -1, -Cxid2 / dxi_f[j-1]);

                //         setCoeff(0, 0, Cxid1 / (2 * dxi_f[j - 1]));
                //         setCoeff(0, -1, -Cxid1 / (2 * dxi_f[j - 1]));
                //     }

                //     continue;
                // } else {
                    // d^2 psi / d p^2
                    // setCoeff(+1, 0, dp2C / dp_f[i]);
                    // setCoeff(-1, 0, dp2C / dp_f[i-1]);
                    // setCoeff( 0, 0, -dp2C * (dp_f[i]+dp_f[i-1]) / (dp_f[i]*dp_f[i-1]));

                    // d^2 psi / d xi^2
                    // setCoeff(0, +1, dxi2C / (dxi_f[j] * dxi_f[j]);
                    // setCoeff(0, -1, dxi2C / dxi_f[j-1]);
                    // setCoeff(0,  0, -dxi2C * (dxi_f[j]+dxi_f[j-1]) / (dxi_f[j]*dxi_f[j-1]));

                    // d psi / d p
                    // setCoeff(+1, 0, dpC / (dp_f[i] + dp_f[i-1]));
                    // setCoeff(-1, 0,-dpC / (dp_f[i] + dp_f[i-1]));

                    // d psi / d xi
                    // setCoeff(0, +1, Cxid1 / (dxi_f[j] + dxi_f[j-1]));
                    // setCoeff(0, -1,-Cxid1 / (dxi_f[j] + dxi_f[j-1]));

                    // // psi
                    // setCoeff(0, 0, 1 - a*a);
                // }
                // if (i == 0) {
                //     setCoeff(0, 0, -3.0/2.0*Cpd1/dp[0]);
                //     setCoeff(+1, 0, 2*Cpd1/dp[0]);
                //     setCoeff(+2, 0, -1.0/2.0*Cpd1/dp[0]);
                //     setCoeff(0, 0, 2*Cpd2/pow(dp[0], 2));
                //     setCoeff(+1, 0, -5*Cpd2/pow(dp[0], 2));
                //     setCoeff(+2, 0, 4*Cpd2/pow(dp[0], 2));
                //     setCoeff(+3, 0, -Cpd2/pow(dp[0], 2));
                // } else {
                //     setCoeff(-1, 0, -1.0/2.0*Cpd1/dp[0]);
                //     setCoeff(0, 0, 0);
                //     setCoeff(1, 0, (1.0/2.0)*Cpd1/dp[0]);
                //     setCoeff(-1, 0, Cpd2/pow(dp[0], 2));
                //     setCoeff(0, 0, -2*Cpd2/pow(dp[0], 2));
                //     setCoeff(1, 0, Cpd2/pow(dp[0], 2));
                // }


                // if (i == 0) {
                //     setCoeff(0, 0, -Cpd1*(dp[i + 1] + 2*dp[i])/(dp[i]*(dp[i + 1] + dp[i])));
                //     setCoeff(+1, 0, Cpd1/dp[i] + Cpd1/dp[i + 1]);
                //     setCoeff(+2, 0, -Cpd1*dp[i]/(dp[i + 1]*(dp[i + 1] + dp[i])));
                //     setCoeff(0, 0, 2*Cpd2*(2*dp[i + 1] + dp[i + 2] + 3*dp[i])/(dp[i]*(pow(dp[i + 1], 2) + dp[i + 1]*dp[i + 2] + 2*dp[i + 1]*dp[i] + dp[i + 2]*dp[i] + pow(dp[i], 2))));
                //     setCoeff(+1, 0, -2*Cpd2*(2*dp[i + 1] + dp[i + 2] + 2*dp[i])/(dp[i + 1]*dp[i]*(dp[i + 1] + dp[i + 2])));
                //     setCoeff(+2, 0, 2*Cpd2*(dp[i + 1] + dp[i + 2] + 2*dp[i])/(dp[i + 1]*dp[i + 2]*(dp[i + 1] + dp[i])));
                //     setCoeff(+3, 0, -2*Cpd2*(dp[i + 1] + 2*dp[i])/(dp[i + 2]*(pow(dp[i + 1], 2) + 2*dp[i + 1]*dp[i + 2] + dp[i + 1]*dp[i] + pow(dp[i + 2], 2) + dp[i + 2]*dp[i])));
                // } else if (i == np1 - 1) {
                //     std::abort();
                //     // setCoeff(-2, 0, Cpd1*dp[i-1]/(dp[i-2]*(dp[i-1] + dp[i-2])));
                //     // setCoeff(-1, 0, -Cpd1/dp[i-2] - Cpd1/dp[i-1]);
                //     // setCoeff(0, 0, Cpd1*(2*dp[i-1] + dp[i-2])/(dp[i-1]*(dp[i-1] + dp[i-2])));
                //     // setCoeff(-3, 0, -2*Cpd2*(2*dp[i-1] + dp[i-2])/(dp[i-3]*(dp[i-1]*dp[i-2] + dp[i-1]*dp[i-3] + pow(dp[i-2], 2) + 2*dp[i-2]*dp[i-3] + pow(dp[i-3], 2))));
                //     // setCoeff(-2, 0, 2*Cpd2*(2*dp[i-1] + dp[i-2] + dp[i-3])/(dp[i-2]*dp[i-3]*(dp[i-1] + dp[i-2])));
                //     // setCoeff(-1, 0, -2*Cpd2*(2*dp[i-1] + 2*dp[i-2] + dp[i-3])/(dp[i-1]*dp[i-2]*(dp[i-2] + dp[i-3])));
                //     // setCoeff(0, 0, 2*Cpd2*(3*dp[i-1] + 2*dp[i-2] + dp[i-3])/(dp[i-1]*(pow(dp[i-1], 2) + 2*dp[i-1]*dp[i-2] + dp[i-1]*dp[i-3] + pow(dp[i-2], 2) + dp[i-2]*dp[i-3])));
                // } else {
                //     setCoeff(-1, 0, -Cpd1*dp[i]/(dp[i-1]*(dp[i-1] + dp[i])));
                //     setCoeff(0, 0, -Cpd1/dp[i] + Cpd1/dp[i-1]);
                //     setCoeff(1, 0, Cpd1*dp[i-1]/(dp[i]*(dp[i-1] + dp[i])));
                //     setCoeff(-1, 0, 2*Cpd2/(dp[i-1]*(dp[i-1] + dp[i])));
                //     setCoeff(0, 0, -2*Cpd2/(dp[i-1]*dp[i]));
                //     setCoeff(1, 0, 2*Cpd2/(dp[i]*(dp[i-1] + dp[i])));
                // }

                // if (j == 0) {
                //     setCoeff(0, +1, Cxid1*dxi0/(dxi[j]*(-dxi0 - dxi[j])));
                //     setCoeff(0, 0, -Cxid1/dxi[j] + Cxid1/dxi0);
                //     setCoeff(0, 0, -Cxid1*dxi[j]/(dxi0*(-dxi0 - dxi[j])));
                //     setCoeff(0, +1, -2*Cxid2/(dxi[j]*(-dxi0 - dxi[j])));
                //     setCoeff(0, 0, 2*Cxid2/(-dxi0*dxi[j]));
                //     setCoeff(0, 0, 2*Cxid2/(-dxi0*(-dxi0 - dxi[j])));
                // } else if (j == np2 - 1) {
                //     setCoeff(0, -1, -Cxid1*dximax/(dxi[j-1]*(dxi[j-1] + dximax)));
                //     setCoeff(0, 0, -Cxid1/dximax + Cxid1/dxi[j-1]);
                //     setCoeff(0, 0, Cxid1*dxi[j-1]/(dximax*(dxi[j-1] + dximax)));
                //     setCoeff(0, -1, 2*Cxid2/(dxi[j-1]*(dxi[j-1] + dximax)));
                //     setCoeff(0, 0, -2*Cxid2/(dxi[j-1]*dximax));
                //     setCoeff(0, 0, 2*Cxid2/(dximax*(dxi[j-1] + dximax)));
                // } else {
                //     setCoeff(0, -1, -Cxid1*dxi[j]/(dxi[j-1]*(dxi[j-1] + dxi[j])));
                //     setCoeff(0, 0, -Cxid1/dxi[j] + Cxid1/dxi[j-1]);
                //     setCoeff(0, +1, Cxid1*dxi[j-1]/(dxi[j]*(dxi[j-1] + dxi[j])));
                //     setCoeff(0, -1, 2*Cxid2/(dxi[j-1]*(dxi[j-1] + dxi[j])));
                //     setCoeff(0, 0, -2*Cxid2/(dxi[j-1]*dxi[j]));
                //     setCoeff(0, +1, 2*Cxid2/(dxi[j]*(dxi[j-1] + dxi[j])));
                // // }

                // setCoeff(0, 0, 1 - a*a);




                const real_t dp2C = gamma[i] * gamma[i] / dp[i];
                const real_t dxi2C = (1 - xi[j] * xi[j]) / (p[i]*p[i] * dxi[j]);
                const real_t dpC = (2 / p[i] + 3*p[i]);
                const real_t dxiC =-2*xi[j]/(p[i]*p[i]);

                if (/*i == 0 ||*/ j == 0 || j == np2 - 1) {
                    /*if (i == 0) {
                        // setCoeff(1, 0, dp2C / dp_f[0]);
                        // setCoeff(0, 0, -dp2C / dp_f[0]);

                        // setCoeff(1, 0, dpC / dp_f[0]);
                        // setCoeff(0, 0, -dpC / dp_f[0]);
                        } else*/ if (j == 0) {
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
                    if (i == 0) {
                        // d^2 psi / d p^2
                        setCoeff( 0, 0, 2 * dp2C / (dp_f[i]));
                        setCoeff( 1, 0, -5 * dp2C / (dp_f[i]));
                        setCoeff( 2, 0, 4 * dp2C / (dp_f[i]));
                        setCoeff( 3, 0, -dp2C / (dp_f[i]));

                        // d psi / d p
                        setCoeff(0, 0, -3.0 / 2.0 * dpC / dp_f[i]);
                        setCoeff(1, 0, 2.0 * dpC / dp_f[i]);
                        setCoeff(2, 0, -1.0/2.0 * dpC / dp_f[i]);
                    } else {
                        // d^2 psi / d p^2
                        setCoeff(+1, 0, dp2C / dp_f[i]);
                        setCoeff(-1, 0, dp2C / dp_f[i-1]);
                        setCoeff( 0, 0, -dp2C * (dp_f[i]+dp_f[i-1]) / (dp_f[i]*dp_f[i-1]));

                        // d psi / d p
                        setCoeff(+1, 0, dpC / (dp_f[i] + dp_f[i-1]));
                        setCoeff(-1, 0,-dpC / (dp_f[i] + dp_f[i-1]));
                    }

                    // d^2 psi / d xi^2
                    setCoeff(0, +1, dxi2C / dxi_f[j]);
                    setCoeff(0, -1, dxi2C / dxi_f[j-1]);
                    setCoeff(0,  0, -dxi2C * (dxi_f[j]+dxi_f[j-1]) / (dxi_f[j]*dxi_f[j-1]));


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
