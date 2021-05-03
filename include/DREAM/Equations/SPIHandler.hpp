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
        real_t dt;

        len_t spi_velocity_mode;
        len_t spi_ablation_mode;
        len_t spi_deposition_mode;
        len_t spi_heat_absorbtion_mode;
        len_t spi_cloud_radius_mode;
        
        real_t VpVolNormFactor;
        real_t rclPrescribedConstant;

        real_t *gradRCartesian=nullptr;
        real_t *gradRCartesianPrevious=nullptr;
        real_t rCoordPClosestApproach;
        len_t nSplit=1;
        len_t iSplit=0;
        bool turningPointPassed;

        len_t id_ncold;
        len_t id_ni;
        len_t id_Tcold;
        len_t id_rp;
        len_t id_xp;
        len_t id_vp;
        len_t id_Wcold;
        len_t id_Whot;
        len_t id_qhot;
        len_t id_ntot;

        real_t *ncold;
        real_t *Tcold;
        real_t *rp;
        real_t *rp_initial;
        real_t *rpPrevious;
        real_t *xpPrevious;
        real_t *xp;
        real_t *vp;
        real_t *Wcold;
        real_t *Whot;
        real_t *qhot;
        real_t *ntot;

        real_t *Ypdot=nullptr;
        real_t *rCld=nullptr;
        real_t *depositionRate=nullptr;
        real_t *depositionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionRate=nullptr;
        real_t *rCoordPPrevious=nullptr;
        real_t *thetaCoordPPrevious=nullptr;
        real_t *rCoordPNext=nullptr;
        real_t *thetaCoordPNext=nullptr;
        len_t *irp=nullptr;
        real_t *qtot=nullptr;
        real_t *Eeff=nullptr;
        real_t *pelletMolarMass=nullptr;
        real_t *pelletMolarVolume=nullptr;
        real_t *pelletDensity=nullptr;
        real_t *lambda=nullptr;
        real_t *NGSConstantFactor=nullptr;

        static const len_t nMolarMassList;
        static const len_t ZMolarMassList[];
        static const len_t isotopesMolarMassList[];
        static const real_t molarMassList[];
        static const len_t nSolidDensityList;
        static const len_t ZSolidDensityList[];
        static const len_t isotopesSolidDensityList[];
        static const real_t solidDensityList[];

        void CalculateYpdotNGSParksTSDW();
        void CalculateYpdotNGSParksTSDWKinetic();
        void CalculateAdiabaticHeatAbsorbtionRateMaxwellian();

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

        void Rebuild(real_t dt);

        void evaluatePartialContributionYpdotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        void evaluatePartialContributionYpdotNGSKinetic(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        void evaluatePartialContributionDepositionRateDensCons(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset);
        void evaluatePartialContributionAdiabaticHeatAbsorbtionRateMaxwellian(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        void evaluatePartialContributionYpdot(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        void evaluatePartialContributionDepositionRate(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset);
        void evaluatePartialContributionAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        real_t *GetYpdot() {return this->Ypdot;}
        real_t *CalculateDepositionRate(real_t *SPIMolarFraction);
        real_t *GetHeatAbsorbtionRate() {return this->heatAbsorbtionRate;}
        
        len_t GetNShard(){return this->nShard;}


    };
}
#endif
