#ifndef _DREAM_EQUATIONS_ION_HANDLER_HPP
#define _DREAM_EQUATIONS_ION_HANDLER_HPP

namespace DREAM { class IonHandler; }

#include "FVM/config.h"
#include "FVM/Grid/RadialGrid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

namespace DREAM {
    class IonHandler{

    private:
        len_t nr,  // number of radial grid points
              nZ,  // number of atomic species
              nzs; // number of ions (including charge states)
        len_t niID;
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        len_t *Zs;
        //len_t *ZList, *Z0List;
        len_t *ZOffsets;
        

        virtual void Initialize();
        virtual void DeallocateAll();
//        virtual void GetOffset(len_t Z, len_t Z0, real_t *&offsets, len_t *nOffsets);
    public:

        IonHandler(FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, len_t *Z, len_t NZ);
        virtual ~IonHandler();
        
        const len_t GetNZ() const { return nZ; }
        const len_t GetNzs() const { return nzs; }

        const len_t* GetZs() const{return Zs;}
        
        const len_t GetIndex(len_t iz, len_t Z0) const{return ZOffsets[iz]+Z0;}

        virtual const real_t GetIonDensityAtZ(len_t ir, len_t Z, len_t Z0) const;
        virtual const real_t GetIonDensity(len_t ir, len_t iz, len_t Z0) const;
        virtual const real_t* GetIonDensity(len_t ir, len_t iZ) const;
        virtual const real_t GetTotalIonDensity(len_t ir, len_t iZ) const;
        virtual const real_t GetTritiumDensity(len_t ir, len_t *tritiumIndices, len_t numTritiumIndices) const;
        virtual const real_t GetFreePlusBoundElectronDensity(len_t ir) const;
        virtual const real_t GetZeff(len_t ir) const;
        virtual const real_t GetZtot(len_t ir) const;
        

    };
}



#endif/*_DREAM_EQUATIONS_ION_HANDLER_HPP*/
