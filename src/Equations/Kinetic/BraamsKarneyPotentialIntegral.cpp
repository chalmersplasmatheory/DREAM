#include "DREAM/Equations/Kinetic/BraamsKarneyPotentialIntegral.hpp"

#include <gsl/gsl_sf_ellint.h>

using namespace DREAM;

// TODO: Move to a better place.
static bool DoIntegralBC(len_t i, len_t j, len_t np1, len_t np2) {
    (void)j, (void)np2;
    return i == np1 - 1;
}

bool BraamsKarneyPotentialIntegral::GridRebuilt() {
    cache.clear();
    return EquationTerm::GridRebuilt();
}

bool BraamsKarneyPotentialIntegral::SetJacobianBlock(const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t *) {
    if (derivId == unknId) {
        this->SetMatrixElements(jac, nullptr);
        return true;
    }

    return false;
}

BraamsKarneyPotentialIntegral::BraamsKarneyPotentialIntegral(FVM::Grid *grid,
                                                             BraamsKarneyPotentialIntegral::Type type)
    : EquationTerm(grid), type(type) {
    SetName("Braams-Karney Potential Integral");

    GridRebuilt();
}

static real_t Gamma(const real_t p) {
    return sqrt(1 + p * p);
}

template<typename F>
static real_t IntegrateFPhi(F f, real_t p, real_t xi, real_t pprime, real_t xiprime) {
    const int N = 24;

    real_t total = 0;

    for (int i = 0; i < N; i++) {
        real_t theta = M_PI * i / (double)(N - 1);
        real_t r = Gamma(p) * Gamma(pprime) - p * pprime * xi * xiprime
            - p * pprime * sqrt((1 - xi * xi) * (1 - xiprime * xiprime)) * cos(theta);

        if (r < 1)
            r = 1;

        // Trapezoidal rule.
        if (i == 0 || i == N - 1)
            total += f(r) / 2;
        else
            total += f(r);
    }

    return 2 * M_PI * total / (N - 1);
}

static real_t EvaluateSpecialFunction(
    const real_t p, const real_t xi,
    const real_t pprime, const real_t xiprime,
    bool calcH
    ) {
    real_t a2 = sqrt((1+p*p)*(1+pprime*pprime)) - p*pprime*xi*xiprime;
    real_t b4 = p*p*pprime*pprime*(1-xi*xi)*(1-xiprime*xiprime);

    // Equivalent to |xi|, |xiprime| = 1
    if (b4 == 0) {
        if (calcH)
            return (2*M_PI*a2 / sqrt(a2*a2 - 1));
        else
            return (2*M_PI    / sqrt(a2*a2 - 1));
    }

    real_t b2 = sqrt(b4);
    real_t a4m1mb4 =
        p*p + pprime*pprime + p*p*pprime*pprime*(xi*xi + xiprime*xiprime) -
        2*sqrt((1+p*p)*(1+pprime*pprime))*p*pprime*xi*xiprime;

    real_t sqr2 = a4m1mb4*a4m1mb4 - 4*b4;
    // Avoid round-off error (if the expression is negative,
    // it is a machine precision error ==> I ~ H ~ 0)
    if (sqr2 <= 0)
        return 0;

    real_t sqr = sqrt(sqr2);
    // lambda+ & lambda-
    real_t lp = (a4m1mb4 + sqr) / (2*b4);
    real_t lm = (a4m1mb4 - sqr) / (2*b4);

    // t+ & t-
    real_t tp = lp*b2*a2 / (1 + lp*b4);
    real_t tm = lm*b2*a2 / (1 + lm*b4);

    // Evaluate elliptic integrals
    real_t k = sqrt(lm/lp);
    real_t K = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);

    // I
    real_t a2p_num = b4*(2 + a4m1mb4 + sqr);
    real_t a2p2_dn = 4*sqr2;
    real_t a1m_num = lp * b4 * (2 + a4m1mb4 - sqr);
    real_t I = 4*K / ((tp - tm)*sqrt(a2p_num*a1m_num / a2p2_dn));

    if (!calcH)
        return I;

    real_t P = gsl_sf_ellint_Pcomp(k, -tm/tp, GSL_PREC_DOUBLE);

    // H
    // a2-b2*tp = a2*lm/(lm+1)
    return (a2*lm/(lm+1)*I + 4*b2 / sqrt(a2p_num*a1m_num / a2p2_dn) * P);
}

real_t BraamsKarneyPotentialIntegral::Integrand(real_t p, real_t xi, real_t pprime, real_t xiprime) {
    switch (type) {
    case BraamsKarneyPotentialIntegral::Type::UPS0:
        return -EvaluateSpecialFunction(p, xi, pprime, xiprime, false) / Gamma(pprime) / (4 * M_PI);
    case BraamsKarneyPotentialIntegral::Type::UPS1:
        return -IntegrateFPhi([] (real_t r) { return sqrt(r * r - 1); }, p, xi, pprime, xiprime)
            / Gamma(pprime) / (8 * M_PI);
    case BraamsKarneyPotentialIntegral::Type::UPS2:
        return -IntegrateFPhi([] (real_t r) { return r * acosh(r) - sqrt(r * r - 1); }, p, xi, pprime, xiprime)
            / Gamma(pprime) / (32 * M_PI);
    case BraamsKarneyPotentialIntegral::Type::PI0:
        return -EvaluateSpecialFunction(p, xi, pprime, xiprime, true) / Gamma(pprime) / (4 * M_PI);
    case BraamsKarneyPotentialIntegral::Type::PI1:
        return -IntegrateFPhi([] (real_t r) { return acosh(r); }, p, xi, pprime, xiprime)
            / Gamma(pprime) / (8 * M_PI);
    default: std::abort();
    }
}

template<typename F1, typename F2>
void BraamsKarneyPotentialIntegral::SetElements(F1 X, F2 ApplyX) {
    const len_t nr  = grid->GetNr();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const real_t
            *dp1 = mg->GetDp1(),
            *dp2 = mg->GetDp2();
        const real_t
            *p1 = mg->GetP1(),
            *p2 = mg->GetP2();

        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            for (len_t j = 0; j < np2; j++) {
                len_t idx = offset + j * np1 + i;

                if (!DoIntegralBC(i, j, np1, np2))
                    continue;

                auto [cached_integrand, has_value] = getCachedIntegrand(ir, i, j);

                if (!has_value) {
                    cached_integrand.resize(np1 * np2);
                    for (len_t ip = 0; ip < np1; ip++){
                        for (len_t jp = 0; jp < np2; jp++){
                            cached_integrand[jp * np1 + ip] =
                                Integrand(p1[i], p2[j], p1[ip], p2[jp]) * p1[ip] * p1[ip] * dp2[jp] * dp1[ip];
                        }
                    }
                }
                for (len_t ip = 0; ip < np1; ip++){
                    for (len_t jp = 0; jp < np2; jp++){
                        len_t idx_p = offset + jp * np1 + ip;
                        X(idx, idx_p, cached_integrand[idx_p - offset]);
                    }
                }

                ApplyX(idx);
            }
        }
        offset += np1*np2;
    }
}

void BraamsKarneyPotentialIntegral::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    len_t N_IND = 0;
    const len_t SIZE = GetNumberOfNonZerosPerRow();
    PetscInt *IND_ARR = new PetscInt[SIZE];
    PetscScalar *VAL_ARR = new PetscScalar[SIZE];
    SetElements([&](auto idx, auto idx_p, auto V) {
                    (void)idx;
                    IND_ARR[N_IND] = idx_p;
                    VAL_ARR[N_IND]= (V);
                    N_IND++;
                },
        [&](auto idx) {
            mat->SetRow(idx, N_IND, IND_ARR, VAL_ARR);
            N_IND = 0;
        }
        );
    delete [] IND_ARR;
    delete [] VAL_ARR;
}

void BraamsKarneyPotentialIntegral::SetVectorElements(real_t *vec, const real_t *f) {
    SetElements([&](auto idx, auto idx_p, auto V) {
                    vec[idx] += f[idx_p] * V;
                },
        [&](auto) {
        }
        );
}
