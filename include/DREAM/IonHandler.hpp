#ifndef _DREAM_EQUATIONS_ION_HANDLER_HPP
#define _DREAM_EQUATIONS_ION_HANDLER_HPP

namespace DREAM { class IonHandler; }

#include <string>
#include <vector>
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
        const len_t *Zs;  // List of atomic charges for each species (size nZ)
        len_t *ZOffsets;

        std::vector<std::string> ionNames;
        
        virtual void DeallocateAll();

        
    public:

        IonHandler(FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, const len_t *Z, len_t NZ, std::vector<std::string>& ionNames);
        virtual ~IonHandler();

        virtual void Initialize(); // Call it rebuild?

        const len_t GetNZ() const { return nZ; }
        const len_t GetNzs() const { return nzs; }

        const len_t  GetZ(const len_t iZ) const { return Zs[iZ]; }
        const len_t* GetZs() const {return Zs;}
        
        const len_t GetIndex(len_t iz, len_t Z0) const{return ZOffsets[iz]+Z0;}

        const std::string& GetName(const len_t iZ) { return this->ionNames[iZ]; }
        const std::vector<std::string>& GetNameList() { return this->ionNames; }

        virtual const real_t GetIonDensityAtZ(len_t ir, len_t Z, len_t Z0) const;
        virtual const real_t GetIonDensity(len_t ir, len_t iz, len_t Z0) const;
        virtual const real_t* GetIonDensity(len_t ir, len_t iZ) const;
        virtual const real_t GetTotalIonDensity(len_t ir, len_t iZ) const;
        virtual const real_t GetTritiumDensity(len_t ir, len_t *tritiumIndices, len_t numTritiumIndices) const;
        virtual real_t* evaluateFreePlusBoundElectronDensityFromQuasiNeutrality(real_t *ntot=nullptr);
        virtual real_t* evaluateZeff();
        virtual real_t* evaluateZtot();
        
    };
}



#endif/*_DREAM_EQUATIONS_ION_HANDLER_HPP*/
