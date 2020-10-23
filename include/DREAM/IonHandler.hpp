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
        std::vector<std::string> tritiumNames;
        len_t nTritium=0;
        int_t *tritiumIndices;
        
        virtual void DeallocateAll();

        
    public:

        IonHandler(FVM::RadialGrid *rg, FVM::UnknownQuantityHandler *u, const len_t *Z, len_t NZ, std::vector<std::string>& ionNames, std::vector<std::string>& tritiumName);
        virtual ~IonHandler();

        virtual void Initialize(); // Call it rebuild?

        const len_t GetNZ() const { return nZ; }
        const len_t GetNzs() const { return nzs; }

        const len_t  GetZ(const len_t iZ) const { return Zs[iZ]; }
        const len_t* GetZs() const {return Zs;}
        
        const len_t GetIndex(len_t iz, len_t Z0) const{return ZOffsets[iz]+Z0;}
        void GetIonIndices(len_t, len_t&, len_t&);


        const std::string& GetName(const len_t iZ) { return this->ionNames[iZ]; }
        const std::vector<std::string>& GetNameList() { return this->ionNames; }
        const std::vector<std::string>& GetTritiumNameList() { return this->tritiumNames; }
        const int_t *GetTritiumIndices() const { return this->tritiumIndices; }
        const len_t GetNTritiumIndices() const { return this->nTritium; }

        const real_t GetIonDensityAtZ(len_t ir, len_t Z, len_t Z0) const;
        const real_t GetIonDensity(len_t ir, len_t iz, len_t Z0) const;
        const real_t* GetIonDensity(len_t ir, len_t iZ) const;
        const real_t GetTotalIonDensity(len_t ir, len_t iZ) const;
        const real_t GetTritiumDensity(len_t ir) const;
        real_t* evaluateFreePlusBoundElectronDensityFromQuasiNeutrality(real_t *ntot=nullptr);
        real_t* evaluateFreeElectronDensityFromQuasiNeutrality(real_t *nfree=nullptr);
        real_t evaluateFreeElectronDensityFromQuasiNeutrality(len_t ir);
        real_t* evaluateBoundElectronDensityFromQuasiNeutrality(real_t *nbound=nullptr);
        real_t evaluateBoundElectronDensityFromQuasiNeutrality(len_t ir);

        real_t* evaluateZeff();
        real_t evaluateZeff(len_t);
        real_t evaluateZ0Z0(len_t);
        real_t *evaluateZeff0();
        real_t evaluateZeff0(len_t);
        real_t *evaluateZ0_Z();
        real_t evaluateZ0_Z(len_t);
        real_t *evaluateZ0Z();
        real_t evaluateZ0Z(len_t);

        real_t* evaluateZtot();
        
    };
}



#endif/*_DREAM_EQUATIONS_ION_HANDLER_HPP*/
