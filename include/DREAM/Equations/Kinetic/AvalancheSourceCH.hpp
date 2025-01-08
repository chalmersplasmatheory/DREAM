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

    private:
        len_t id_fgrid, id_fre;
        len_t id_Efield;
        real_t scaleFactor;
        real_t preFactor;
        real_t pCutoff;
        CHSourceMode sourceMode;
        bool isRunawayGrid;
        FVM::Grid* runawayGrid;
        
        int_t *pIn_Indices;
        std::vector<len_t> *f_pIndices;
        std::vector<len_t> *f_xiIndices;
        std::vector<len_t> *f_re_pIndices;
        std::vector<len_t> *f_re_xiIndices;

        // If the quantity is a hot-tail source, but there is a runaway 
        // grid availabel, this will be true
        bool htgridWithREgrid = false; 

        void SetPinIndices();

    public:
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;

        AvalancheSourceCH(FVM::Grid*, FVM::UnknownQuantityHandler*, real_t, real_t, CHSourceMode sm = CH_SOURCE_MODE_KINETIC, bool isRunawayGrid=true, FVM::Grid* runawayGrid=nullptr);
        ~AvalancheSourceCH();

        real_t EvaluateCHSource(len_t ir, len_t i, len_t j);
    };
}

#endif/*_DREAM_EQUATIONS_AVALANCHE_SOURCE_CH_HPP*/


