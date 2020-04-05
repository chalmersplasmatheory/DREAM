#ifndef _DREAM_CONSTANTS_HPP
#define _DREAM_CONSTANTS_HPP

#include "FVM/config.h"

namespace DREAM {
    class Constants{

    public:

                
        static const real_t pi; 

        //Speed of light in m/s:
        static const real_t c; 

        //Classical electron radius in m:
        static const real_t r0; 

        //Elementary charge in C:
        static const real_t ec;

        //Electron mass in kg:
        static const real_t me;

        //Electron rest energy in eV:
        static const real_t mc2inEV;

        //Deuterium mass in kg
        static const real_t mD; 

        //Vacuum permittivity in SI:
        static const real_t eps0;

        //Vacuum permeability in SI:
        static const real_t mu0;

        //Planck's reduced constant in Js
        static const real_t hbar;

        //Compton wavelength in m
        static const real_t lambda_c;

        //Bohr radius in m
        static const real_t a0;                

        //Fine structure constant
        static const real_t alpha;
        

        

    };
}

#endif/*_DREAM_CONSTANTS_HPP*/
