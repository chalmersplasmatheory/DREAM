#ifndef _DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP
#define _DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include <limits>
namespace DREAM {
    class AvalancheSourceCH
        : public FluidSourceTerm {
    public:
        enum CHSourceMode{
            CH_SOURCE_MODE_FLUID,
            CH_SOURCE_MODE_KINETIC
        };

        // Setting to determine whether to apply the source term for
        // negative or positive xi. In the "adaptive" case, this is determined
        // by the sign of the electric field.
        /*enum CHSourcePitchMode {
            CH_SOURCE_PITCH_ADAPTIVE,
            CH_SOURCE_PITCH_POSITIVE,
            CH_SOURCE_PITCH_NEGATIVE
        };*/
    private:
        //len_t id_ntot;
        len_t id_fgrid, id_fre;
        len_t id_Efield;
        real_t scaleFactor;
        real_t preFactor;
        real_t pCutoff;
        //CHSourcePitchMode sourceXiMode;
        CHSourceMode sourceMode;
        bool isRunawayGrid;
        FVM::Grid* runawayGrid;
        
        int_t *pIn_Indices;
        std::vector<len_t> *f_pIndices;
        std::vector<len_t> *f_xiIndices;
        std::vector<len_t> *f_re_pIndices;
        std::vector<len_t> *f_re_xiIndices;

        bool htgridWithREgrid = false; // If the quantity is a hot-tail source, but there is a runaway grid availabel, this will be true

        void SetPinIndices();

    //protected: //TODO: ok?
    public:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    //public: // TODO: above
        AvalancheSourceCH(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t, real_t, /*CHSourcePitchMode sxm = CH_SOURCE_PITCH_ADAPTIVE,*/ CHSourceMode sm = CH_SOURCE_MODE_KINETIC, bool isRunawayGrid=true, FVM::Grid* runawayGrid=nullptr);
        ~AvalancheSourceCH();

        real_t EvaluateCHSource(len_t ir, len_t i, len_t j);
    };
}

#endif/*_DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP*/


