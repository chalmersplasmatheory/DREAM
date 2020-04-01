#ifndef _DREAM_CONSTANTS_HPP
#define _DREAM_CONSTANTS_HPP

#include "FVM/config.h"

namespace DREAM {
    class Constants{

    public:
        static const real_t pi;
    
        //Speed of light in m/s:
        static const real_t c;
    
        //Vacuum permittivity in SI:
        static const real_t eps0;
    
        //Electron mass in kg:
        static const real_t m_e;

        //Classical electron radius in m:
        static const real_t r0;

        //Elementary charge in C:
        static const real_t ec;

        

        

    };
}

#endif/*_DREAM_CONSTANTS_HPP*/
