/**
 * Header file for a class that calculates and stores quantities related to the SPI shards
 */
#ifndef _DREAM_EQUATIONS_SPI_HANDLER_HPP
#define _DREAM_EQUATIONS_SPI_HANDLER_HPP

namespace DREAM { class SPIHandler; }
#include <iostream>

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM{
    class SPIHandler{
    private:
        //FVM::ScalarGrid *sGrid;
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        len_t nShard;
        len_t nr;

        len_t spi_velocity_mode;
        len_t spi_ablation_mode;
        len_t spi_deposition_mode;
        len_t spi_heat_absorbtion_mode;
        len_t spi_cloud_radius_mode;

        real_t pelletMolarMass;
        real_t pelletMolarVolume;
        real_t pelletDensity;
        real_t lambda;
        real_t VpVolNormFactor;
        real_t rclPrescribedConstant;

        len_t id_ncold;
        len_t id_ntot;
        len_t id_ni;
        len_t id_Tcold;
        len_t id_rp;
        len_t id_xp;
        len_t id_vp;

        real_t *ncold;
        real_t *ntot;
        real_t *Tcold;
        real_t *rp;
        real_t *rp_initial;
        //real_t *rpPrevious;
        real_t *xpPrevious;
        real_t *xp;
        real_t *vp;

        real_t *rpdot=nullptr;
        real_t *rCld=nullptr;
        real_t *depositionRate=nullptr;
        real_t *depositionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionRate=nullptr;
        real_t *rCoordPPrevious=nullptr;
        real_t *rCoordPNext=nullptr;
        len_t *irp=nullptr;

        static const len_t nMolarMassList;
        static const len_t ZMolarMassList[];
        static const len_t isotopesMolarMassList[];
        static const real_t molarMassList[];
        static const len_t nSolidDensityList;
        static const len_t ZSolidDensityList[];
        static const len_t isotopesSolidDensityList[];
        static const real_t solidDensityList[];

        void CalculateRpDotNGSParksTSDW();
        void CalculateDepositionRate();
        void CalculateAdiabaticHeatAbsorbtionRateNGS();

        real_t rSourceMax; 
        real_t rSourceMin;
        void CalculateTimeAveragedDeltaSourceLocal(real_t *timeAveragedDeltaSource);
        void CalculateGaussianSourceLocal(real_t *gaussianSource);
        real_t CalculateRDotDepositionLocal(len_t ir);

        void CalculateIrp();
        void CalculateRCld();
        real_t CalculateLambda(real_t X);

    public:
        SPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t *Z, len_t *isotopes, const real_t *molarFraction, len_t NZ, 
            OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode,
            OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode,
            OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode,
            OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode,
            OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode, real_t VpVolNormFactor, real_t rclPrescribedConstant);
        ~SPIHandler();
        void AllocateQuantities();
        void DeallocateQuantities();

        void Rebuild();

        void evaluatePartialContributionRpDotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        void evaluatePartialContributionDepositionRateNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor, real_t SPIMolarFraction, len_t rOffset);
        void evaluatePartialContributionAdiabaticHeatAbsorbtionRateNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        void evaluatePartialContributionRpDot(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        void evaluatePartialContributionDepositionRate(FVM::Matrix *jac,len_t derivId, real_t scaleFactor, real_t SPIMolarFraction, len_t rOffset);
        void evaluatePartialContributionAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        real_t *GetRpdot() {return this->rpdot;}
        real_t *GetDepositionRate() {
            return this->depositionRate;}
        real_t *GetHeatAbsorbtionRate() {return this->heatAbsorbtionRate;}


    };
}
#endif
