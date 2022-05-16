#include "DREAM/Equations/Kinetic/BraamsKarneyDiffusion.hpp"

using namespace DREAM;

void BraamsKarneyDiffusion::DCoeffs::init(len_t nr, len_t *n1, len_t *n2) {
    this->dd11 = new real_t*[nr];
    this->dd12 = new real_t*[nr];
    this->dd21 = new real_t*[nr];
    this->dd22 = new real_t*[nr];

    len_t
        nElements_f1 = 0,
        nElements_f2 = 0;

    for (len_t i = 0; i < nr; i++) {
        nElements_f1 += (n1[i]+1)*n2[i];
        nElements_f2 += n1[i]*(n2[i]+1);
    }

	this->dd11[0] = new real_t[nElements_f1];
	this->dd12[0] = new real_t[nElements_f1];
	this->dd22[0] = new real_t[nElements_f2];
	this->dd21[0] = new real_t[nElements_f2];

	for (len_t i = 1; i < nr; i++) {
		this->dd11[i] = this->dd11[i-1] + ((n1[i-1]+1)*n2[i-1]);
		this->dd12[i] = this->dd12[i-1] + ((n1[i-1]+1)*n2[i-1]);
		this->dd22[i] = this->dd22[i-1] + (n1[i-1]*(n2[i-1]+1));
		this->dd21[i] = this->dd21[i-1] + (n1[i-1]*(n2[i-1]+1));
	}
}

BraamsKarneyDiffusion::DCoeffs::~DCoeffs() {
	if (dd11 != nullptr) {
        delete [] dd11[0];
        delete [] dd11;
    }
    if (dd12 != nullptr) {
        delete [] dd12[0];
        delete [] dd12;
    }
    if (dd21 != nullptr) {
        delete [] dd21[0];
        delete [] dd21;
    }
    if (dd22 != nullptr) {
        delete [] dd22[0];
        delete [] dd22;
    }
}

/**
 * Constructor.
 */
BraamsKarneyDiffusion::BraamsKarneyDiffusion(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, CoulombLogarithm *lnLambda, len_t id_upsilon_1, len_t id_upsilon_2)
	: FVM::DiffusionTerm(g, unknowns), lnLambda(lnLambda), id_upsilon_1(id_upsilon_1), id_upsilon_2(id_upsilon_2) {

    SetName("BraamsKarneyDiffusion");

    AddUnknownForJacobian(unknowns, id_upsilon_1);
    AddUnknownForJacobian(unknowns, id_upsilon_2);

	for (int i = 0; i < dds_n1; i++)
		for (int j = 0; j < dds_n2; j++)
		dds[j][i].init(nr, n1, n2);
}

const real_t alphaBarOverLnLambda = 3.75927427447396469133e-19; // ec^4 / (eps0^2 * me^2 * c^3)

template<typename F>
void setFD_XI(F setCoeff, int i, int j, int np1, int np2, const real_t *dp, const real_t *dxi, real_t dximax, real_t dxi0, real_t diff_p1, real_t diff_xi1p1, real_t diff_xi1, real_t diff_xi2) {
// d/dp
    if (diff_p1 != 0.0) {
        if (i == 0) {
            if (j == 0) {
                setCoeff(0, 0, -3.0/4.0*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(1, 0, diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(2, 0, -1.0/4.0*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(0, 0, (3.0/4.0)*diff_p1*dxi[0]/(dp[0]*dxi0));
                setCoeff(1, 0, -diff_p1*dxi[0]/(dp[0]*dxi0));
                setCoeff(2, 0, (1.0/4.0)*diff_p1*dxi[0]/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(0, -1, -9.0/4.0*diff_p1/dp[0]);
                setCoeff(1, -1, 3*diff_p1/dp[0]);
                setCoeff(2, -1, -3.0/4.0*diff_p1/dp[0]);
                setCoeff(0, -1, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(1, -1, -diff_p1/dp[0]);
                setCoeff(2, -1, (1.0/4.0)*diff_p1/dp[0]);
            } else {
                setCoeff(0, 0, -9.0/4.0*diff_p1/dp[0]);
                setCoeff(1, 0, 3*diff_p1/dp[0]);
                setCoeff(2, 0, -3.0/4.0*diff_p1/dp[0]);
                setCoeff(0, -1, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(1, -1, -diff_p1/dp[0]);
                setCoeff(2, -1, (1.0/4.0)*diff_p1/dp[0]);
            }
        } else if (i == np1 - 1) {
            if (j == 0) {
                setCoeff(-2, 0, (1.0/4.0)*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(-1, 0, -diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(0, 0, (3.0/4.0)*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(-2, 0, -1.0/4.0*diff_p1*dxi[0]/(dp[0]*dxi0));
                setCoeff(-1, 0, diff_p1*dxi[0]/(dp[0]*dxi0));
                setCoeff(0, 0, -3.0/4.0*diff_p1*dxi[0]/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(-2, -1, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(-1, -1, -3*diff_p1/dp[0]);
                setCoeff(0, -1, (9.0/4.0)*diff_p1/dp[0]);
                setCoeff(-2, -1, -1.0/4.0*diff_p1/dp[0]);
                setCoeff(-1, -1, diff_p1/dp[0]);
                setCoeff(0, -1, -3.0/4.0*diff_p1/dp[0]);
            } else {
                setCoeff(-2, 0, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(-1, 0, -3*diff_p1/dp[0]);
                setCoeff(0, 0, (9.0/4.0)*diff_p1/dp[0]);
                setCoeff(-2, -1, -1.0/4.0*diff_p1/dp[0]);
                setCoeff(-1, -1, diff_p1/dp[0]);
                setCoeff(0, -1, -3.0/4.0*diff_p1/dp[0]);
            }
        } else {
            if (j == 0) {
                setCoeff(-1, 0, -1.0/4.0*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(1, 0, (1.0/4.0)*diff_p1*(2*dxi0 + dxi[0])/(dp[0]*dxi0));
                setCoeff(-1, 0, (1.0/4.0)*diff_p1*dxi[0]/(dp[0]*dxi0));
                setCoeff(1, 0, -1.0/4.0*diff_p1*dxi[0]/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(-1, -1, -3.0/4.0*diff_p1/dp[0]);
                setCoeff(1, -1, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(-1, -1, (1.0/4.0)*diff_p1/dp[0]);
                setCoeff(1, -1, -1.0/4.0*diff_p1/dp[0]);
            } else {
                setCoeff(-1, 0, -3.0/4.0*diff_p1/dp[0]);
                setCoeff(1, 0, (3.0/4.0)*diff_p1/dp[0]);
                setCoeff(-1, -1, (1.0/4.0)*diff_p1/dp[0]);
                setCoeff(1, -1, -1.0/4.0*diff_p1/dp[0]);
            }
        }
    }
// d2/dxidp
    if (diff_xi1p1 != 0.0) {
        if (i == 0) {
            if (j == 0) {
                setCoeff(0, 0, -3.0/2.0*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(1, 0, 2*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(2, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(0, 0, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(1, 0, -2*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(2, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(0, -1, -3.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, 2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(2, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, -2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(2, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
            } else {
                setCoeff(0, 0, -3.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, 0, 2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(2, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, -2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(2, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
            }
        } else if (i == np1 - 1) {
            if (j == 0) {
                setCoeff(-2, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(-1, 0, -2*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(0, 0, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(-2, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(-1, 0, 2*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(0, 0, -3.0/2.0*diff_xi1p1/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(-2, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, -2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-2, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, 2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, -3.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            } else {
                setCoeff(-2, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, 0, -2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, 0, (3.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-2, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, 2*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, -3.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            }
        } else {
            if (j == 0) {
                setCoeff(-1, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(1, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(-1, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi0));
                setCoeff(1, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi0));
            } else if (j == np2 - 0) {
                setCoeff(-1, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            } else {
                setCoeff(-1, 0, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, 0, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(1, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            }
        }
    }
// d2/dxi2
    if (diff_xi2 != 0.0) {
        if (j == 0) {
            setCoeff(0, 1, diff_xi2*(4*dxi0 + 5*dxi[0])/(dxi[0]*(pow(dxi0, 2) + 3*dxi0*dxi[0] + 2*pow(dxi[0], 2))));
            setCoeff(0, 0, -diff_xi2*(4*dxi0 + 3*dxi[0])/(dxi0*dxi[0]*(dxi0 + dxi[0])));
            setCoeff(0, 0, diff_xi2*(2*dxi0 + 3*dxi[0])/(dxi0*dxi[0]*(dxi0 + dxi[0])));
            setCoeff(0, 1, -diff_xi2*(2*dxi0 + dxi[0])/(dxi[0]*(pow(dxi0, 2) + 3*dxi0*dxi[0] + 2*pow(dxi[0], 2))));
        } else if (j == 1) {
            setCoeff(0, 1, (1.0/2.0)*diff_xi2*(2*dxi0 + 7*dxi[0])/(pow(dxi[0], 2)*(dxi0 + 2*dxi[0])));
            setCoeff(0, 0, -diff_xi2*(2*dxi0 + 5*dxi[0])/(pow(dxi[0], 2)*(dxi0 + dxi[0])));
            setCoeff(0, -1, diff_xi2/pow(dxi[0], 2) + (3.0/2.0)*diff_xi2/(dxi0*dxi[0]));
            setCoeff(0, -1, -3*diff_xi2*dxi[0]/(dxi0*(pow(dxi0, 2) + 3*dxi0*dxi[0] + 2*pow(dxi[0], 2))));
        } else if (j == np2 - 0) {
            setCoeff(0, -2, (3.0/2.0)*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -1, -7.0/2.0*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -1, (5.0/2.0)*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -2, -1.0/2.0*diff_xi2/pow(dxi[0], 2));
        } else if (j == np2 - 1) {
            setCoeff(0, 0, 9*diff_xi2*dxi[0]/(dximax*(2*pow(dxi[0], 2) + 3*dxi[0]*dximax + pow(dximax, 2))));
            setCoeff(0, 0, -9.0/2.0*diff_xi2/(dxi[0]*dximax) + diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -1, diff_xi2*(7*dxi[0] - 2*dximax)/(pow(dxi[0], 2)*(dxi[0] + dximax)));
            setCoeff(0, -2, (1.0/2.0)*diff_xi2*(-5*dxi[0] + 2*dximax)/(pow(dxi[0], 2)*(2*dxi[0] + dximax)));
        } else {
            setCoeff(0, 1, (3.0/2.0)*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, 0, -7.0/2.0*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -1, (5.0/2.0)*diff_xi2/pow(dxi[0], 2));
            setCoeff(0, -2, -1.0/2.0*diff_xi2/pow(dxi[0], 2));
        }
    }
// d/dxi
    if (diff_xi1 != 0.0) {
        if (j == 0) {
            setCoeff(0, 1, diff_xi1/dxi[0]);
            setCoeff(0, 0, -diff_xi1/dxi[0]);
        } else if (j == 1) {
            setCoeff(0, 1, diff_xi1/dxi[0]);
            setCoeff(0, 0, -diff_xi1/dxi[0]);
        } else if (j == np2 - 0) {
            setCoeff(0, -2, diff_xi1/dxi[0]);
            setCoeff(0, -1, -diff_xi1/dxi[0]);
        } else if (j == np2 - 1) {
            setCoeff(0, 0, 2*diff_xi1*dxi[0]/(dximax*(dxi[0] + dximax)));
            setCoeff(0, 0, -2*diff_xi1/dximax + diff_xi1/dxi[0]);
            setCoeff(0, -1, diff_xi1*(dxi[0] - dximax)/(dxi[0]*(dxi[0] + dximax)));
        } else {
            setCoeff(0, 1, diff_xi1/dxi[0]);
            setCoeff(0, 0, -diff_xi1/dxi[0]);
        }
    }
}

template<typename F>
void setFD_P(F setCoeff, int i, int j, int np1, int np2, const real_t *dp, const real_t *dxi, real_t dximax, real_t dxi0, real_t diff_p1, real_t diff_p2, real_t diff_xi1p1, real_t diff_xi1) {
// d/dp
    if (diff_p1 != 0.0) {
        if (i == 1) {
            setCoeff(-1, 0, -23.0/24.0*diff_p1/dp[0]);
            setCoeff(0, 0, (7.0/8.0)*diff_p1/dp[0]);
            setCoeff(1, 0, (1.0/8.0)*diff_p1/dp[0]);
            setCoeff(2, 0, -1.0/24.0*diff_p1/dp[0]);
        } else if (i == np1) {
            setCoeff(-4, 0, -23.0/24.0*diff_p1/dp[0]);
            setCoeff(-3, 0, (31.0/8.0)*diff_p1/dp[0]);
            setCoeff(-2, 0, -47.0/8.0*diff_p1/dp[0]);
            setCoeff(-1, 0, (71.0/24.0)*diff_p1/dp[0]);
        } else if (i == np1 - 1) {
            setCoeff(-3, 0, (1.0/24.0)*diff_p1/dp[0]);
            setCoeff(-2, 0, -1.0/8.0*diff_p1/dp[0]);
            setCoeff(-1, 0, -7.0/8.0*diff_p1/dp[0]);
            setCoeff(0, 0, (23.0/24.0)*diff_p1/dp[0]);
        } else {
            setCoeff(-2, 0, (1.0/24.0)*diff_p1/dp[0]);
            setCoeff(-1, 0, -9.0/8.0*diff_p1/dp[0]);
            setCoeff(0, 0, (9.0/8.0)*diff_p1/dp[0]);
            setCoeff(1, 0, -1.0/24.0*diff_p1/dp[0]);
        }
    }
// d2/dp2
    if (diff_p2 != 0.0) {
        if (i == 1) {
            setCoeff(-1, 0, (3.0/2.0)*diff_p2/pow(dp[0], 2));
            setCoeff(0, 0, -7.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(1, 0, (5.0/2.0)*diff_p2/pow(dp[0], 2));
            setCoeff(2, 0, -1.0/2.0*diff_p2/pow(dp[0], 2));
        } else if (i == np1) {
            setCoeff(-4, 0, -3.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(-3, 0, (11.0/2.0)*diff_p2/pow(dp[0], 2));
            setCoeff(-2, 0, -13.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(-1, 0, (5.0/2.0)*diff_p2/pow(dp[0], 2));
        } else if (i == np1 - 1) {
            setCoeff(-3, 0, -1.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(-2, 0, (5.0/2.0)*diff_p2/pow(dp[0], 2));
            setCoeff(-1, 0, -7.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(0, 0, (3.0/2.0)*diff_p2/pow(dp[0], 2));
        } else {
            setCoeff(-2, 0, (1.0/2.0)*diff_p2/pow(dp[0], 2));
            setCoeff(-1, 0, -1.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(0, 0, -1.0/2.0*diff_p2/pow(dp[0], 2));
            setCoeff(1, 0, (1.0/2.0)*diff_p2/pow(dp[0], 2));
        }
    }
// d/dxi
    if (diff_xi1 != 0.0) {
        if (i == np1) {
            if (j == 0) {
                setCoeff(-2, 1, -1.0/2.0*diff_xi1*dxi0/(dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-1, 1, (3.0/2.0)*diff_xi1*dxi0/(dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-2, 0, (1.0/2.0)*diff_xi1*dxi[0]/(dxi0*(dxi0 + dxi[0])));
                setCoeff(-1, 0, -3.0/2.0*diff_xi1*dxi[0]/(dxi0*(dxi0 + dxi[0])));
                setCoeff(-2, 0, (1.0/2.0)*diff_xi1*(dxi0 - dxi[0])/(dxi0*dxi[0]));
                setCoeff(-1, 0, (3.0/2.0)*diff_xi1*(-dxi0 + dxi[0])/(dxi0*dxi[0]));
            } else if (j == np2 - 1) {
                setCoeff(-2, 0, -1.0/2.0*diff_xi1*dxi[0]/(dximax*(dxi[0] + dximax)));
                setCoeff(-1, 0, (3.0/2.0)*diff_xi1*dxi[0]/(dximax*(dxi[0] + dximax)));
                setCoeff(-2, -1, (1.0/2.0)*diff_xi1*dximax/(dxi[0]*(dxi[0] + dximax)));
                setCoeff(-1, -1, -3.0/2.0*diff_xi1*dximax/(dxi[0]*(dxi[0] + dximax)));
                setCoeff(-2, 0, (1.0/2.0)*diff_xi1*(dxi[0] - dximax)/(dxi[0]*dximax));
                setCoeff(-1, 0, (3.0/2.0)*diff_xi1*(-dxi[0] + dximax)/(dxi[0]*dximax));
            } else {
                setCoeff(-2, 1, -1.0/4.0*diff_xi1/dxi[0]);
                setCoeff(-1, 1, (3.0/4.0)*diff_xi1/dxi[0]);
                setCoeff(-2, -1, (1.0/4.0)*diff_xi1/dxi[0]);
                setCoeff(-1, -1, -3.0/4.0*diff_xi1/dxi[0]);
            }
        } else {
            if (j == 0) {
                setCoeff(-1, 1, (1.0/2.0)*diff_xi1*dxi0/(dxi[0]*(dxi0 + dxi[0])));
                setCoeff(0, 1, (1.0/2.0)*diff_xi1*dxi0/(dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-1, 0, -1.0/2.0*diff_xi1*dxi[0]/(dxi0*(dxi0 + dxi[0])));
                setCoeff(0, 0, -1.0/2.0*diff_xi1*dxi[0]/(dxi0*(dxi0 + dxi[0])));
                setCoeff(-1, 0, (1.0/2.0)*diff_xi1*(-dxi0 + dxi[0])/(dxi0*dxi[0]));
                setCoeff(0, 0, (1.0/2.0)*diff_xi1*(-dxi0 + dxi[0])/(dxi0*dxi[0]));
            } else if (j == np2 - 1) {
                setCoeff(-1, 0, (1.0/2.0)*diff_xi1*dxi[0]/(dximax*(dxi[0] + dximax)));
                setCoeff(0, 0, (1.0/2.0)*diff_xi1*dxi[0]/(dximax*(dxi[0] + dximax)));
                setCoeff(-1, -1, -1.0/2.0*diff_xi1*dximax/(dxi[0]*(dxi[0] + dximax)));
                setCoeff(0, -1, -1.0/2.0*diff_xi1*dximax/(dxi[0]*(dxi[0] + dximax)));
                setCoeff(-1, 0, (1.0/2.0)*diff_xi1*(-dxi[0] + dximax)/(dxi[0]*dximax));
                setCoeff(0, 0, (1.0/2.0)*diff_xi1*(-dxi[0] + dximax)/(dxi[0]*dximax));
            } else {
                setCoeff(-1, 1, (1.0/4.0)*diff_xi1/dxi[0]);
                setCoeff(0, 1, (1.0/4.0)*diff_xi1/dxi[0]);
                setCoeff(-1, -1, -1.0/4.0*diff_xi1/dxi[0]);
                setCoeff(0, -1, -1.0/4.0*diff_xi1/dxi[0]);
            }
        }
    }
// d2/dxidp
    if (diff_xi1p1 != 0.0) {
        if (i == np1) {
            if (j == 0) {
                setCoeff(-2, 1, -diff_xi1p1*dxi0/(dp[0]*dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-1, 1, diff_xi1p1*dxi0/(dp[0]*dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-2, 0, diff_xi1p1*dxi[0]/(dp[0]*dxi0*(dxi0 + dxi[0])));
                setCoeff(-1, 0, -diff_xi1p1*dxi[0]/(dp[0]*dxi0*(dxi0 + dxi[0])));
                setCoeff(-2, 0, diff_xi1p1*(dxi0 - dxi[0])/(dp[0]*dxi0*dxi[0]));
                setCoeff(-1, 0, diff_xi1p1*(-dxi0 + dxi[0])/(dp[0]*dxi0*dxi[0]));
            } else if (j == np2 - 1) {
                setCoeff(-2, 0, -diff_xi1p1*dxi[0]/(dp[0]*dximax*(dxi[0] + dximax)));
                setCoeff(-1, 0, diff_xi1p1*dxi[0]/(dp[0]*dximax*(dxi[0] + dximax)));
                setCoeff(-2, -1, diff_xi1p1*dximax/(dp[0]*dxi[0]*(dxi[0] + dximax)));
                setCoeff(-1, -1, -diff_xi1p1*dximax/(dp[0]*dxi[0]*(dxi[0] + dximax)));
                setCoeff(-2, 0, diff_xi1p1*(dxi[0] - dximax)/(dp[0]*dxi[0]*dximax));
                setCoeff(-1, 0, diff_xi1p1*(-dxi[0] + dximax)/(dp[0]*dxi[0]*dximax));
            } else {
                setCoeff(-2, 1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, 1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-2, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            }
        } else {
            if (j == 0) {
                setCoeff(-1, 1, -diff_xi1p1*dxi0/(dp[0]*dxi[0]*(dxi0 + dxi[0])));
                setCoeff(0, 1, diff_xi1p1*dxi0/(dp[0]*dxi[0]*(dxi0 + dxi[0])));
                setCoeff(-1, 0, diff_xi1p1*dxi[0]/(dp[0]*dxi0*(dxi0 + dxi[0])));
                setCoeff(0, 0, -diff_xi1p1*dxi[0]/(dp[0]*dxi0*(dxi0 + dxi[0])));
                setCoeff(-1, 0, diff_xi1p1*(dxi0 - dxi[0])/(dp[0]*dxi0*dxi[0]));
                setCoeff(0, 0, diff_xi1p1*(-dxi0 + dxi[0])/(dp[0]*dxi0*dxi[0]));
            } else if (j == np2 - 1) {
                setCoeff(-1, 0, -diff_xi1p1*dxi[0]/(dp[0]*dximax*(dxi[0] + dximax)));
                setCoeff(0, 0, diff_xi1p1*dxi[0]/(dp[0]*dximax*(dxi[0] + dximax)));
                setCoeff(-1, -1, diff_xi1p1*dximax/(dp[0]*dxi[0]*(dxi[0] + dximax)));
                setCoeff(0, -1, -diff_xi1p1*dximax/(dp[0]*dxi[0]*(dxi[0] + dximax)));
                setCoeff(-1, 0, diff_xi1p1*(dxi[0] - dximax)/(dp[0]*dxi[0]*dximax));
                setCoeff(0, 0, diff_xi1p1*(-dxi[0] + dximax)/(dp[0]*dxi[0]*dximax));
            } else {
                setCoeff(-1, 1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, 1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(-1, -1, (1.0/2.0)*diff_xi1p1/(dp[0]*dxi[0]));
                setCoeff(0, -1, -1.0/2.0*diff_xi1p1/(dp[0]*dxi[0]));
            }
        }
    }

}

template<typename T1, typename T2>
void BraamsKarneyDiffusion::SetCoefficients(T1 psi, T2 phi, bool overwrite,
											real_t **d11, real_t **d12,
											real_t **d21, real_t **d22) {
	len_t r_offset = 0;
	const len_t nr = grid->GetNr();

	for (len_t ir = 0; ir < nr; ir++) {
		auto *mg = grid->GetMomentumGrid(ir);
		const len_t np1 = mg->GetNp1(), np2 = mg->GetNp2();        

        CoulombLogarithm::collqty_settings settings;
        settings.lnL_type = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
        auto alphabar = [&] (real_t p) {
                            return alphaBarOverLnLambda * lnLambda->evaluateAtP(ir, p, &settings);
                            // return alphaBarOverLnLambda * lnLambda->evaluateLnLambdaT(ir);
                            // return alphaBarOverLnLambda * lnLambda->evaluateLnLambdaT(ir);
                        };

		const real_t *p_g = mg->GetP1(),
			*p_f = mg->GetP1_f(),
			*xi_f = mg->GetP2_f(),
			*dp_f = mg->GetDp1_f(),
			*xi_g = mg->GetP2(),
			*dp = mg->GetDp1(),
			*dxi = mg->GetDp2();

        real_t dxi0 = 2 * (xi_g[0] + 1);
        real_t dximax = 2 * (1 - xi_g[np2 - 1]);

        auto callPhi = [&](len_t i, len_t j, int di, int dj) {
                           len_t idx = r_offset + (i + di) + np1 * (j + dj);
                           return phi(idx, di , dj);
                       };

        auto callPsi = [&](len_t i, len_t j, int di, int dj) {
                           len_t idx = r_offset + (i + di) + np1 * (j + dj);
                           return psi(idx, di , dj);
                       };

        for (len_t j = 0; j < np2; j++) 
			for (len_t i = 1; i < np1 + 1; i++) {
                real_t d12_val = 0;
                real_t d11_val = 0;
                real_t p = i == np1 ? p_f[i - 1] + dp_f[i - 1] : p_f[i];
                //real_t sintheta = sqrt(1 - xi_g[j] * xi_g[j]);
                real_t gamma = sqrt(1 + p * p);
                setFD_P([&] (int di, int dj, real_t coeff) {
                            d11_val += coeff * callPsi(i, j, di, dj);
                        }, i, j, np1, np2, dp, dxi, dximax, dxi0,
                    p, gamma * gamma, 0, 0
                    );
                d11_val *= gamma;

                setFD_P([&] (int di, int dj, real_t coeff) {
                            d12_val += coeff * callPsi(i, j, di, dj);
                        }, i, j, np1, np2, dp, dxi, dximax, dxi0,
                    0, 0, 1, -1 / p
                    );
                d12_val *= gamma * (1 - xi_g[j] * xi_g[j]) / (p * p);

                if (i == np1) {
                    d11_val += gamma * callPhi(i, j, -1, 0);
                } else {
                    d11_val += gamma * (callPhi(i, j, 0, 0) +
                                                callPhi(i, j, -1, 0)) / 2;
                }

                d11_val *= -alphabar(p_f[i]);
                d12_val *= -alphabar(p_f[i]);

                if (overwrite) {
                    D11(ir, i, j, d11) = d11_val;
                    D12(ir, i, j, d12) = d12_val;
                } else {
                    D11(ir, i, j, d11) += d11_val;
                    D12(ir, i, j, d12) += d12_val;
                }
            }

        for (len_t j = 0; j < np2 + 1; j++) 
			for (len_t i = 0; i < np1; i++) {
                real_t d22_val = 0;
                real_t d21_val = 0;
                real_t p = p_g[i];
                real_t xi = j == np2 ? 1 : xi_f[j];
                //real_t sintheta = sqrt(1 - xi * xi);
                real_t gamma = sqrt(1 + p * p);
                setFD_XI([&] (int di, int dj, real_t coeff) {
                            d21_val += coeff * callPsi(i, j, di, dj);
                        }, i, j, np1, np2, dp, dxi, dximax, dxi0,
                    0, 1, -1 / p, 0
                    );

                d21_val *= gamma * (1 - xi * xi) / (p * p);

                setFD_XI([&] (int di, int dj, real_t coeff) {
                            d22_val += coeff * callPsi(i, j, di, dj);
                        }, i, j, np1, np2, dp, dxi, dximax, dxi0,
                    1 / (p * p) + 1, 0, -xi / (p * p * p), (1 - xi * xi) / (p * p * p)
                    );

                d22_val *= (1 - xi * xi) / (p * gamma);

                real_t plus_norm = (1 - xi * xi) / (gamma * p * p);
                if (j == np2) {
                    d22_val += plus_norm * callPhi(i, j, 0, -1);
                } else if (j == 0) {
                    d22_val += plus_norm * callPhi(i, j, 0, 0);
                } else {
                    d22_val += plus_norm * (callPhi(i, j, 0, 0) +
                                            callPhi(i, j, 0, -1)) / 2;
                }

                d21_val *= -alphabar(p_g[i]);
                d22_val *= -alphabar(p_g[i]);

                if (overwrite) {
                    D21(ir, i, j, d21) = d21_val;
                    D22(ir, i, j, d22) = d22_val;
                } else {
                    D21(ir, i, j, d21) += d21_val;
                    D22(ir, i, j, d22) += d22_val;
                }
            }

		r_offset += np1 * np2;
	}
}

void BraamsKarneyDiffusion::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler *x
){
    real_t *ups_1 = x->GetUnknownData(id_upsilon_1);
    real_t *ups_2 = x->GetUnknownData(id_upsilon_2);

    SetCoefficients(
        [&] (len_t idx, int, int) {
            return ups_1[idx] - 4 * ups_2[idx];
        },
        [&] (len_t idx, int, int) {
            return ups_1[idx] + 4 * ups_2[idx];
        }, false, this->d11, this->d12, this->d21, this->d22);
}

// Set jacobian of the diffusion coefficients for this diffusion term
void BraamsKarneyDiffusion::SetPartialDiffusionTerm(len_t derivId, len_t) {
	for (int i = 0; i < dds_n1; i++) {
		for (int j = 0; j < dds_n2; j++) {
			real_t psi_coeff = 1, phi_coeff = 1;
			if (derivId == id_upsilon_2) {
				psi_coeff = -4;
				phi_coeff = 4;
			} else if (derivId != id_upsilon_1) {
				std::abort();
			}
			SetCoefficients(
				[&] (len_t, int di, int dj) {
					return (di == n1_offsets[i] && dj == n2_offsets[j]) ? psi_coeff : 0;
				},
				[&] (len_t, int di, int dj) {
					return (di == n1_offsets[i] && dj == n2_offsets[j]) ? phi_coeff : 0;
				}, true, this->dds[j][i].dd11, this->dds[j][i].dd12, this->dds[j][i].dd21, this->dds[j][i].dd22);
		}
	}
}

// Overridden
void BraamsKarneyDiffusion::SetPartialJacobianContribution(int_t diagonalOffset, jacobian_interp_mode, len_t, FVM::Matrix *jac, const real_t *x, bool){
	if (diagonalOffset)
		return;

	for (int k1 = 0; k1 < dds_n1; k1++) {
		for (int k2 = 0; k2 < dds_n2; k2++) {
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
#   include "BraamsKarneyDiffusion.set.cpp"
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
