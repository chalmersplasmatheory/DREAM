#include<iostream>
#include<cmath>
#include "DREAM/Constants.hpp"

using namespace std;


real_t aiCoefficient(int_t i, int_t Z, int_t N) {
    real_t c[4][5] = {
        {1.1831, 0.1738, 0.0913, 0.0182, 0.7702},
        {0.8368, 1.0987, 0.9642, 1.2535, 0.2618},
        {0.3841, 0.6170, 1.0000, 1.0000, 1.0000},
        {0.5883, 0.0461, 1.0000, 1.0000, 1.0000}
    };
    real_t lambda = c[0][i] * pow(Z, c[1][i]);
    real_t n = c[2][i] * pow(Z, c[3][i]);
    real_t x = (real_t)(Z - N)/ (real_t)Z;
    real_t ai = 1.0 / sqrt(lambda * lambda  * ( 1.0 - pow(x, n + 1)) / (1.0 - x));
    return ai;
}

real_t yiCoefficient(int_t i, int_t Z, int_t N, real_t p) {
    real_t ai = aiCoefficient(i, Z, N);
    real_t yi = 2 * ai * p / DREAM::Constants::alpha;
    return yi;
}

real_t gEquation(int_t Z, int_t Z0, real_t p) {
    int_t N = Z-Z0;

    int_t Ni[5] = { 
        std::min(N, (int_t)2),
        std::min( std::max(N-2, (int_t)0) , (int_t)8),
        std::min( std::max(N-10, (int_t)0) , (int_t)18),
        std::min( std::max(N-28, (int_t)0), (int_t)26),
        std::max(N-54, (int_t)0)
    };

    int_t Zi[5];
    for(int i=0;i<5;i++) {
        Zi[i] = Z - Ni[i];
    }

    real_t y[5];
    for(int i=0;i<5;i++) {
        real_t yi = yiCoefficient(i, Z, N, p);
        y[i] = yi * yi;
    }

    real_t a[5];
    for(int i=0;i<5;i++) {
        real_t ai = aiCoefficient(i, Z, N);
        a[i] = ai * ai;
    }

    real_t lnyi[5];
    for(int i=0;i<5;i++) {
        lnyi[i] = log(1 + y[i]);
    }

    real_t sum = 0.0;
    for(int i=0;i<5;i++) {
        sum+=((Z*Z - Zi[i]*Zi[i])* lnyi[i] - Ni[i] * Ni[i] *  y[i] / ( 1.0 + y[i]))/2;
    }

    for(int i=0;i<5;i++)
		for(int k=i+1;k<5;k++) {
			sum-= Ni[i] * Ni[k] * (lnyi[i] + a[k] / (a[k] - a[i]) * (lnyi[k] - lnyi[i]));
		}
    return sum;
}

