#include "DREAM/Constants.hpp"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
using namespace DREAM;

//Speed of light in m/s:
const real_t Constants::c = 299792458.0; 

//Classical electron radius in m:
const real_t Constants::r0 = 2.8179403227e-15; 

//Elementary charge in C:
const real_t Constants::ec = 1.60217662e-19;

//Electron mass in kg:
const real_t Constants::me = 9.10938356e-31;

//Electron rest energy in eV:
const real_t Constants::mc2inEV = 0.51099895000e6;

//Deuterium mass in kg
const real_t Constants::mD = 3.3435837724e-27; 

//Vacuum permittivity in SI:
const real_t Constants::eps0 = 8.85418782e-12;

//Vacuum permeability in SI:
const real_t Constants::mu0 = 1.25663706e-6;

//Planck's reduced constant in Js
const real_t Constants::hbar = 1.054571817e-34;

//Compton wavelength in m
const real_t Constants::lambda_c = 2.4263102367e-12;

//Bohr radius in m
const real_t Constants::a0 = 5.29177210903e-11;

//Fine structure constant
const real_t Constants::alpha = 0.0072973525693;

//Rydberg constant in eV 
const real_t Constants::Ry = 13.605693123; 


//Relativistic Maxwell-Jüttner distribution function at momentum p, density n, temperature T.
//dFdn and dFdT contains the derivative of the distribution function wrt n and T, respectively.
const real_t Constants::RelativisticMaxwellian(const real_t p, const real_t n, const real_t T, real_t *dFdn, real_t *dFdT){
        real_t Theta  = T / mc2inEV;
        real_t K2scaled = gsl_sf_bessel_Knu_scaled(2.0, 1.0/Theta);
        real_t tK2exp = 4*M_PI*Theta * K2scaled;

        const real_t g = sqrt(1+p*p);
        const real_t gMinus1 = p*p/(g+1); // = g-1, for numerical stability for arbitrarily small p
        const real_t e = exp(-gMinus1/Theta);
        
        // set density derivative if requested
        if(dFdn != nullptr)
                *dFdn = 1.0/ tK2exp * e;
        // set temperature derivative if requested
        if(dFdT != nullptr){
                const real_t dedT = gMinus1/(Theta*T)*e;   
                real_t h = Theta * 1e-6;
                real_t K2scaledH = gsl_sf_bessel_Knu_scaled(2.0, 1.0/(Theta+h));
                real_t dtK2dT = 4*M_PI/mc2inEV *( K2scaled + Theta * (K2scaledH-K2scaled)/h);
                *dFdT = n / tK2exp * dedT - n / (tK2exp*tK2exp) * dtK2dT * e;
        }

        return n / tK2exp * e;
}
