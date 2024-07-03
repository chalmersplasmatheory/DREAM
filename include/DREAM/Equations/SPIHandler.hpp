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
#include "DREAM/DREAMException.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM{
    class SPIHandler{
    private:
        //FVM::ScalarGrid *sGrid;
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        RunawayFluid *rf = nullptr;
        len_t nShard;
        len_t nr;
        real_t dt;
        len_t NZ;

        len_t spi_velocity_mode;
        len_t spi_ablation_mode;
        len_t spi_deposition_mode;
        len_t spi_heat_absorbtion_mode;
        len_t spi_cloud_radius_mode;
        len_t spi_magnetic_field_dependence_mode;
        len_t spi_shift_mode;
        
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
        len_t id_Yp;
        len_t id_xp;
        len_t id_vp;
        len_t id_Wcold;
        len_t id_Whot;
        len_t id_qhot;
        len_t id_ntot;

        real_t *ncold;
        real_t *Tcold;
        real_t *Yp;
        real_t *Yp_initial;
        real_t *YpPrevious;
        real_t *xpPrevious;
        real_t *xp;
        real_t *vp;
        real_t *Wcold;
        real_t *Whot;
        real_t *qhot;
        real_t *ntot;
        real_t *ni;


        real_t *Ypdot=nullptr;
        real_t *rCld=nullptr;
        real_t *depositionRate=nullptr;
        real_t *depositionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionProfilesAllShards=nullptr;
        real_t *heatAbsorbtionRate=nullptr;
        real_t *rCoordPPrevious=nullptr;
        real_t *thetaCoordPPrevious=nullptr;
        real_t *phiCoordPPrevious=nullptr;
        real_t *rCoordPNext=nullptr;
        real_t *thetaCoordPNext=nullptr;
        real_t *phiCoordPNext=nullptr;
        len_t *irp=nullptr;
        real_t *qtot=nullptr;
        real_t *Eeff=nullptr;
        real_t *pelletMolarMass=nullptr;
        real_t *pelletMolarVolume=nullptr;
        real_t *pelletDensity=nullptr;
        real_t *lambda=nullptr;
        real_t *NGSConstantFactor=nullptr;
        int_t *nbrShiftGridCell=nullptr;
        int_t *nbrShiftGridCellPrescribed=nullptr;
        real_t *shift_store=nullptr;
        real_t *ncoldPrevious=nullptr;
        real_t *TcoldPrevious=nullptr;
        real_t *YpdotPrevious=nullptr;
        real_t *plasmoidAbsorbtionFactor=nullptr;
        real_t *cosThetaDrift=nullptr;

        static const len_t nMolarMassList;
        static const len_t ZMolarMassList[];
        static const len_t isotopesMolarMassList[];
        static const real_t molarMassList[];
        static const len_t nSolidDensityList;
        static const len_t ZSolidDensityList[];
        static const len_t isotopesSolidDensityList[];
        static const real_t solidDensityList[];

        //Parameters for the shift calculation
        real_t *TDrift=nullptr;
        real_t *ZavgDriftArray;
        real_t *pelletDeuteriumFraction=nullptr;
        real_t* rp=nullptr;
        real_t* rpdot=nullptr;
        real_t* shift_r=nullptr;
        real_t* ZavgDrift=nullptr;
        len_t* isotopesDrift=nullptr;
        real_t T0Drift, DeltaYDrift, RmDrift;
        real_t tAccDrift, tPolDrift, tPeDrift, tExpDrift, tPolDriftPrime, tPeDriftPrime, tExpDriftPrime;
        real_t v0Drift, neBgDrift, niBgDrift, TeBgDrift, BBgDrift, sigmaBgDrift, qBgDrift;
        real_t CSTDrift, CST0Drift, G, n0Drift, a0Drift, tDetachDrift, LcInitDrift, nBarDrift, vLabInitDrift, ReffBgDrift;
        
        // We store the time at the end of the time step for the previous iteration 
        // to be able to tell if the current iteration is the first one on the current time step
        // This is useful for e.g. the drift calculation, which is only rebuilt in the first iteration
        // of every time step to avoid a discontinuity in the jacobian when using a delta fucntion kernel
        real_t t_old=0;

        void CalculateYpdotNGSParksTSDW();
        void CalculateYpdotNGSParksTSDWKinetic();
        real_t CalculateBFieldDampingJOREK(len_t ir);
        void CalculateAdiabaticHeatAbsorbtionRateMaxwellian();

        real_t rSourceMax; 
        real_t rSourceMin;
        void CalculateTimeAveragedDeltaSourceLocal(real_t *timeAveragedDeltaSource);
        void CalculateGaussianSourceLocal(real_t *gaussianSource);
        real_t CalculateRDotDepositionLocal(len_t ir);

        void CalculateIrp();
        int_t CalculateDriftIrp(len_t ip, real_t shift);
        void CalculateRCld();
        real_t CalculateLambda(real_t X);

    public:
        SPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t *Z, len_t *isotopes, const real_t *molarFraction, len_t NZ, 
            OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode,
            OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode,
            OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode,
            OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode,
            OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode, 
            OptionConstants::eqterm_spi_magnetic_field_dependence_mode spi_magnetic_field_dependence_mode, 
            OptionConstants::eqterm_spi_shift_mode spi_shift_mode, 
            real_t *TDrift, real_t T0Drift, real_t DeltaYDrift,real_t RmDrift, real_t *ZavgDriftArray, len_t nZavgDrift, len_t *ZsDrift, len_t *isotopesDrift,
            real_t VpVolNormFactor, real_t rclPrescribedConstant, int_t *nbrShiftGridCell);
        ~SPIHandler();
        void AllocateQuantities();
        void DeallocateQuantities();

        struct integrand_struct {real_t a;};
        // Functions for the drift calculation
        void YpConversion(len_t ip);
        void AssignShardSpecificDriftParameters(len_t ip);
        void AssignDriftComputationParameters(len_t ip);
        void AssignDriftTimeParameters(len_t ip);
        static real_t Integrand(real_t x, void * params);
        real_t Epsiloni(real_t a, real_t b);
        real_t BisFunction(real_t t_prim);
        real_t PrimitiveFirstRow(real_t t_prim);
        real_t PrimitiveSecondRow(real_t t_prim);
        real_t PrimitiveThirdRow(real_t t_prim);
        real_t FirstRow();
        real_t SecondRow();
        real_t ThirdRow();
        real_t Deltar(len_t ip);
        
        void SetREFluid(RunawayFluid *REF) 
        { this->rf = REF; }

        real_t* GetDrift()
        { return this->shift_store; }

        void Rebuild(real_t dt, real_t t);

        bool setJacobianYpdotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        bool setJacobianYpdotNGSKinetic(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        bool setJacobianDepositionRateDensCons(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset);
        bool setJacobianAdiabaticHeatAbsorbtionRateMaxwellian(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        bool setJacobianYpdot(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);
        bool setJacobianDepositionRate(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset);
        bool setJacobianAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac,len_t derivId, real_t scaleFactor);

        real_t *GetYpdot() {return this->Ypdot;}
        real_t *CalculateDepositionRate(real_t *SPIMolarFraction);
        real_t *GetHeatAbsorbtionRate() {return this->heatAbsorbtionRate;}
        
        len_t GetNShard(){return this->nShard;}


    };
}
#endif
