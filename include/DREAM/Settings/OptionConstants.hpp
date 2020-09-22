#ifndef _DREAM_OPTION_CONSTANTS_HPP
#define _DREAM_OPTION_CONSTANTS_HPP

namespace DREAM {
    class OptionConstants {
    public:
        #include "OptionConstants.enum.hpp"

        // CONSTANTS
        // When adding a new constant here, remember
        // to also define it in 'src/Settings/Constants.cpp'.
        // Please, also maintain alphabetical order.
        static const char
            *UQTY_E_FIELD,
            *UQTY_F_HOT,
            *UQTY_F_RE,
            *UQTY_ION_SPECIES,
            *UQTY_I_WALL,
            *UQTY_I_P,
            *UQTY_J_HOT,
            *UQTY_J_OHM,
            *UQTY_J_RE,
            *UQTY_J_TOT,
            *UQTY_N_COLD,
            *UQTY_N_HOT,
            *UQTY_N_RE,
            *UQTY_N_TOT,
            *UQTY_POL_FLUX,
            *UQTY_PSI_WALL,
            *UQTY_PSI_EDGE,
            *UQTY_S_PARTICLE,
            *UQTY_T_COLD,
            *UQTY_V_LOOP_WALL,
            *UQTY_W_COLD;

        // Descriptions of unknown quantities
        static const char
            *UQTY_E_FIELD_DESC,
            *UQTY_F_HOT_DESC,
            *UQTY_F_RE_DESC,
            *UQTY_ION_SPECIES_DESC,
            *UQTY_I_WALL_DESC,
            *UQTY_I_P_DESC,
            *UQTY_J_HOT_DESC,
            *UQTY_J_OHM_DESC,
            *UQTY_J_RE_DESC,
            *UQTY_J_TOT_DESC,
            *UQTY_N_COLD_DESC,
            *UQTY_N_HOT_DESC,
            *UQTY_N_RE_DESC,
            *UQTY_N_TOT_DESC,
            *UQTY_POL_FLUX_DESC,
            *UQTY_PSI_WALL_DESC,
            *UQTY_PSI_EDGE_DESC,
            *UQTY_S_PARTICLE_DESC,
            *UQTY_T_COLD_DESC,
            *UQTY_V_LOOP_WALL_DESC,
            *UQTY_W_COLD_DESC;
    };
}

#endif/*_DREAM_OPTION_CONSTANTS_HPP*/
