#ifndef _DREAM_SETTINGS_EQUATIONS_KINETIC_QUASILINEAR_DIFFUSION_BUILDER_HPP
#define _DREAM_SETTINGS_EQUATIONS_KINETIC_QUASILINEAR_DIFFUSION_BUILDER_HPP

#include <string>
#include "DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"

namespace DREAM {
    class Settings;
    
    /**
     * Builder function for constructing quasilinear diffusion equation terms.
     * 
     * This function encapsulates all the logic for creating QuasilinearDiffusionTerm
     * objects, supporting both pre-computed matrix mode and on-the-fly computation mode.
     * 
     * @param s      Settings object to load configuration from
     * @param mod    Module name (e.g., "eqsys/f_hot")
     * @param grid   Momentum grid on which the distribution function lives
     * @return       Pointer to created QuasilinearDiffusionTerm, or nullptr if disabled
     */
    QuasilinearDiffusionTerm *ConstructQuasilinearDiffusionTerm(
        Settings *s, 
        const std::string& mod, 
        FVM::Grid *grid
    );
}

#endif /* _DREAM_SETTINGS_EQUATIONS_KINETIC_QUASILINEAR_DIFFUSION_BUILDER_HPP */
