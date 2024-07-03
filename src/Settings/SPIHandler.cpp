#include "DREAM/Settings/SimulationGenerator.hpp"

using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/spi"
#define MODULENAME_IONS "eqsys/n_i"

void SimulationGenerator::DefineOptions_SPI(Settings *s){
    s->DefineSetting(MODULENAME "/velocity","method to use for calculating the velocity of the spi shards",(int_t)OptionConstants::EQTERM_SPI_VELOCITY_MODE_NONE);
    s->DefineSetting(MODULENAME "/ablation","method to use for calculating the pellet shard ablation",(int_t)OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/deposition","method to use for calculating the pellet shard deposition",(int_t)OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/heatAbsorbtion","method to use for calculating the heat absorbtion in the neutral cloud surrounding the pellet shards",(int_t)OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/cloudRadiusMode","method to use for calculating the size of the neutral cloud surrounding the pellet shards",(int_t)OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/magneticFieldDependenceMode","method to use for calculating the magnetic field dependence of the ablation rate",(int_t)OptionConstants::EQTERM_SPI_MAGNETIC_FIELD_DEPENDENCE_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/abl_ioniz","method to use for calculating the charge state distribution with which the recently ablated material is deposited",(int_t)OptionConstants::EQTERM_SPI_MAGNETIC_FIELD_DEPENDENCE_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/shift","method to use for displacing the plasma due to the drift",(int_t)OptionConstants::EQTERM_SPI_SHIFT_MODE_NEGLECT);
    s->DefineSetting(MODULENAME "/TDrift","Cloud temperature during the majority of the drift",0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/T0Drift","Cloud temperature directly after the neutral phase", (real_t)0);
    s->DefineSetting(MODULENAME "/DeltaYDrift","Cloud half-width during the drift", (real_t)0);
    s->DefineSetting(MODULENAME "/RmDrift","Major radius to be used for the drift calculation (needed if the grid has an infinite major radius)", (real_t)0);
    s->DefineSetting(MODULENAME "/ZavgDriftArray","Average charge of the ions inside the drifting ablation cloud", 0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/ZsDrift","Ion charges correspnding to ZavgDriftArray", 0, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/isotopesDrift","Ion isotopes correspnding to ZavgDriftArray", 0, (int_t*)nullptr);


    s->DefineSetting(MODULENAME "/init/rp", "initial number of shard particles",0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/init/xp", "initial shard positions",0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/init/vp", "shard velocities",0, (real_t*)nullptr);
    s->DefineSetting(MODULENAME "/init/t_delay", "time delay before the shards start moving",0, (real_t*)nullptr);

    s->DefineSetting(MODULENAME "/VpVolNormFactor", "Norm factor for VpVol=1/R to be used when having an otherwise cylindrical geometry, to get a finita volume of the flux tubes with the correct unit",1.0);
    s->DefineSetting(MODULENAME "/rclPrescribedConstant", "Precribed, constant radius for the neutral cloud surrounding the pellet shards",0.01);
    s->DefineSetting(MODULENAME "/nbrShiftGridCell", "Number of grid cells to shift the deposition", 0, (int_t*)nullptr);
}

SPIHandler *SimulationGenerator::ConstructSPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, Settings *s){
    enum OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode= (enum OptionConstants::eqterm_spi_velocity_mode)s->GetInteger(MODULENAME "/velocity");
    enum OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode = (enum OptionConstants::eqterm_spi_ablation_mode)s->GetInteger(MODULENAME "/ablation");
    enum OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode = (enum OptionConstants::eqterm_spi_deposition_mode)s->GetInteger(MODULENAME "/deposition");
    enum OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode = (enum OptionConstants::eqterm_spi_heat_absorbtion_mode)s->GetInteger(MODULENAME "/heatAbsorbtion");
    enum OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode = (enum OptionConstants::eqterm_spi_cloud_radius_mode)s->GetInteger(MODULENAME "/cloudRadiusMode");
    enum OptionConstants::eqterm_spi_magnetic_field_dependence_mode spi_magnetic_field_dependence_mode = (enum OptionConstants::eqterm_spi_magnetic_field_dependence_mode)s->GetInteger(MODULENAME "/magneticFieldDependenceMode");
    enum OptionConstants::eqterm_spi_shift_mode spi_shift_mode = (enum OptionConstants::eqterm_spi_shift_mode)s->GetInteger(MODULENAME "/shift");

    len_t nZ;
    len_t nZSPInShard;
    len_t nShard;
    len_t nZavgDrift;
    
    const int_t *_Z  = s->GetIntegerArray(MODULENAME_IONS "/Z", 1, &nZ);
    const int_t *_isotopes  = s->GetIntegerArray(MODULENAME_IONS "/isotopes", 1, &nZ);
    const real_t *molarFraction  = s->GetRealArray(MODULENAME_IONS "/SPIMolarFraction", 1, &nZSPInShard);
    real_t VpVolNormFactor = s->GetReal(MODULENAME "/VpVolNormFactor");
    real_t rclPrescribedConstant = s->GetReal(MODULENAME "/rclPrescribedConstant");
    const int_t *_nbrShiftGridCell  = s->GetIntegerArray(MODULENAME "/nbrShiftGridCell", 1, &nShard);
    real_t T0Drift = s->GetReal(MODULENAME "/T0Drift");
    real_t DeltaYDrift = s->GetReal(MODULENAME "/DeltaYDrift");
    real_t RmDrift = s->GetReal(MODULENAME "/RmDrift");

    // Data type conversion
    len_t *Z = new len_t[nZ];
    len_t *isotopes = new len_t[nZ];
    for (len_t i = 0; i < nZ; i++){
        Z[i] = (len_t)_Z[i];
        isotopes[i]=(len_t)_isotopes[i];
    }
    
    int_t *nbrShiftGridCell = new int_t[nShard];
    for (len_t i = 0; i < nShard; i++)
        nbrShiftGridCell[i] = (int_t)_nbrShiftGridCell[i];

    const real_t *_TDrift = s->GetRealArray(MODULENAME "/TDrift", 1, &nShard);
    const real_t *_ZavgDriftArray = s->GetRealArray(MODULENAME "/ZavgDriftArray", 1, &nZavgDrift);
    const int_t *_ZsDrift = s->GetIntegerArray(MODULENAME "/ZsDrift", 1, &nZavgDrift);
    const int_t *_isotopesDrift = s->GetIntegerArray(MODULENAME "/isotopesDrift", 1, &nZavgDrift);
    real_t *TDrift = new real_t[nShard];
    real_t *ZavgDriftArray = new real_t[nZavgDrift];
    len_t *ZsDrift = new len_t[nZavgDrift];
    len_t *isotopesDrift = new len_t[nZavgDrift];
    if(spi_shift_mode == OptionConstants::EQTERM_SPI_SHIFT_MODE_ANALYTICAL){
        for (len_t i = 0; i < nShard; i++)
            TDrift[i] = (real_t)_TDrift[i];
        for (len_t i = 0; i < nZavgDrift; i++){
            ZavgDriftArray[i] = (real_t)_ZavgDriftArray[i];
            ZsDrift[i] = (len_t)_ZsDrift[i];
            isotopesDrift[i] = (len_t)_isotopesDrift[i];
        }
    }

    SPIHandler *SPI=new SPIHandler(g, unknowns, Z, isotopes, molarFraction, nZ, spi_velocity_mode, spi_ablation_mode, spi_deposition_mode, spi_heat_absorbtion_mode, spi_cloud_radius_mode, spi_magnetic_field_dependence_mode, spi_shift_mode, TDrift, T0Drift, DeltaYDrift, RmDrift, ZavgDriftArray, nZavgDrift, ZsDrift, isotopesDrift, VpVolNormFactor, rclPrescribedConstant, nbrShiftGridCell);
    return SPI;
}
