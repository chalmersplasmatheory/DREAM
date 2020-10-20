/**
 * Implementation of a class that calculates and stores quantities related to the SPI shards
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include "DREAM/Equations/SPIHandler.hpp"

using namespace DREAM;
using namespace std;

const len_t SPIHandler::nMolarMassList=3;
const len_t SPIHandler::ZMolarMassList[nMolarMassList]={1,1,10};
const len_t SPIHandler::isotopesMolarMassList[nMolarMassList]={2,0,0};// 0 means naturally occuring mix
const real_t SPIHandler::molarMassList[nMolarMassList]={0.0020141,0.001008,0.020183};// kg/mol

const len_t SPIHandler::nSolidDensityList=3;
const len_t SPIHandler::ZSolidDensityList[nSolidDensityList]={1,1,10};
const len_t SPIHandler::isotopesSolidDensityList[nSolidDensityList]={2,0,0};
const real_t SPIHandler::solidDensityList[nSolidDensityList]={205.9,86,1444};// kg/m^3

const real_t cutoffFrac=0.01;

SPIHandler::SPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t *Z, len_t *isotopes, const real_t *molarFraction, len_t NZ, 
    OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode,
    OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode,
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode,
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode,
    OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode, real_t VpVolNormFactor=1, real_t rclPrescribedConstant=0.01){

    this->rGrid=g->GetRadialGrid();
    //this->sGRid=g->GetScalarGrid();
    this->unknowns=u;
    this->VpVolNormFactor=VpVolNormFactor;

    this->spi_velocity_mode=spi_velocity_mode;
    this->spi_ablation_mode=spi_ablation_mode;
    this->spi_deposition_mode=spi_deposition_mode;
    this->spi_heat_absorbtion_mode=spi_heat_absorbtion_mode;
    this->spi_cloud_radius_mode=spi_cloud_radius_mode;

    if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT){
        this->rclPrescribedConstant=rclPrescribedConstant;
    }else{
        this->rclPrescribedConstant=0.0;
    }

    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_rp = unknowns->GetUnknownID(OptionConstants::UQTY_R_P);
    id_xp = unknowns->GetUnknownID(OptionConstants::UQTY_X_P);
    id_vp = unknowns->GetUnknownID(OptionConstants::UQTY_V_P);


    this->nr=rGrid->GetNr();
    this->nShard=unknowns->GetUnknown(id_rp)->NumberOfMultiples();

    AllocateQuantities();

    real_t molarMass;
    real_t solidDensity;
    real_t pelletDeuteriumFraction=0;
    pelletMolarMass=0;
    pelletMolarVolume=0;
    for(len_t iZ=0;iZ<NZ;iZ++){
        if(molarFraction[iZ]>0){
            for(len_t i=0;i<nMolarMassList;i++){
                if(Z[iZ]==ZMolarMassList[i] && isotopes[iZ]==isotopesMolarMassList[i]){
                    molarMass=molarMassList[i];
                }
            }
            for(len_t i=0;i<nSolidDensityList;i++){
                if(Z[iZ]==ZSolidDensityList[i] && isotopes[iZ]==isotopesSolidDensityList[i]){
                    solidDensity=solidDensityList[i];
                }
            }
            pelletMolarMass+=molarMass*molarFraction[iZ];
            pelletMolarVolume+=molarMass/solidDensity*molarFraction[iZ];
     
            if(Z[iZ]==1 && isotopes[iZ]==2){
                pelletDeuteriumFraction+=molarFraction[iZ];
            }
        }
    }
    pelletDensity=pelletMolarMass/pelletMolarVolume;

    // It seems that the lambda implemented here is only valid for composite neon-deuterium pellets
    // but since the only reference for it is Parks 2017 TSDW presentation it is rather unclear.
    // Also note that lambda in Parks TSDW presentation is defined in terms of the molar fraction of D_2, 
    // while the input gives the molar fraction of D, hence the seemingly weird input argument.
    lambda=CalculateLambda(pelletDeuteriumFraction/2.0/(1.0-pelletDeuteriumFraction/2.0));

    cout<<rGrid->GetVpVol(this->nr-1)<<endl;
}

SPIHandler::~SPIHandler(){
   DeallocateQuantities();
}

void SPIHandler::AllocateQuantities(){
    DeallocateQuantities();

    rpdot = new real_t[nShard];
    rCld = new real_t[nShard];
    depositionRate = new real_t[nr];
    depositionProfilesAllShards = new real_t[nr*nShard];
    heatAbsorbtionRate = new real_t[nr];
    heatAbsorbtionProfilesAllShards = new real_t[nr*nShard];
    rCoordPPrevious = new real_t[nShard];
    rCoordPNext = new real_t[nShard];
    irp = new len_t[nShard];
}

void SPIHandler::DeallocateQuantities(){
    delete [] rpdot;
    delete [] rCld;
    delete [] depositionRate;
    delete [] depositionProfilesAllShards;
    delete [] heatAbsorbtionRate;
    delete [] heatAbsorbtionProfilesAllShards;
    delete [] rCoordPPrevious;
    delete [] rCoordPNext;
    delete [] irp;
}

void SPIHandler::Rebuild(){

    xpPrevious=unknowns->GetUnknownDataPrevious(id_xp);
    xp=unknowns->GetUnknownData(id_xp);
    vp=unknowns->GetUnknownData(id_vp);
    ncold=unknowns->GetUnknownData(id_ncold);
    Tcold=unknowns->GetUnknownData(id_Tcold);
    rp=unknowns->GetUnknownData(id_rp);
    rp_initial=unknowns->GetUnknownInitialData(id_rp);
    //rpPrevious=unknowns->GetUnknownDataPrevious(id_rp);

    if(spi_velocity_mode==OptionConstants::EQTERM_SPI_VELOCITY_MODE_PRESCRIBED){
        // Note that at the first iteration, xp in the current time step will be equal to xpPrevious, unless xp is prescribed!
        for(len_t ip=0;ip<nShard;ip++){
            rCoordPPrevious[ip]=rGrid->GetRFromCartesian(xpPrevious[3*ip],xpPrevious[3*ip+1],xpPrevious[3*ip+2]);
            rCoordPNext[ip]=rGrid->GetRFromCartesian(xp[3*ip],xp[3*ip+1],xp[3*ip+2]);
        }
        CalculateIrp();
    }else if(spi_velocity_mode==OptionConstants::EQTERM_SPI_VELOCITY_MODE_NONE){
        for(len_t ip=0;ip<nShard;ip++){
            rCoordPPrevious[ip]=0;
            rCoordPNext[ip]=0;
            irp[ip]=0;
        }
    }// else {exception}

    
   
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        CalculateRpDotNGSParksTSDW();
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        for(len_t ip=0;ip<nShard;ip++)
            rpdot[ip]=0;
    }// else {exception}

    if(spi_cloud_radius_mode!=OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_NEGLECT)
        CalculateRCld();

    if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL){
        CalculateTimeAveragedDeltaSourceLocal(depositionProfilesAllShards);
        CalculateDepositionRate();
    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){
        CalculateTimeAveragedDeltaSourceLocal(depositionProfilesAllShards);
        for(len_t ip=0;ip<nShard;ip++){
            if(rCoordPNext[ip]>rCoordPPrevious[ip]){
                for(len_t ir=0;ir<nr-1;ir++){
                    depositionProfilesAllShards[ir*nShard+ip]=depositionProfilesAllShards[(ir+1)*nShard+ip];
                }
            }else if(rCoordPNext[ip]<rCoordPPrevious[ip]){
                for(len_t ir=nr-1;ir>0;ir--){
                    depositionProfilesAllShards[ir*nShard+ip]=depositionProfilesAllShards[(ir-1)*nShard+ip];
                }
            }
        }
        CalculateDepositionRate();
    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        CalculateGaussianSourceLocal(depositionProfilesAllShards);
        CalculateDepositionRate();
    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT){
        for(len_t ir=0;ir<nr;ir++)
            depositionRate[ir]=0;
    }// else {exception}

    if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS){
        CalculateTimeAveragedDeltaSourceLocal(heatAbsorbtionProfilesAllShards);
        CalculateAdiabaticHeatAbsorbtionRateNGS();
    }else if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN){
        CalculateGaussianSourceLocal(heatAbsorbtionProfilesAllShards);
        CalculateAdiabaticHeatAbsorbtionRateNGS();
    }else if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT){
        for(len_t ir=0;ir<nr;ir++)
            heatAbsorbtionRate[ir]=0;
    }// else {exception}

}

void SPIHandler::CalculateRpDotNGSParksTSDW(){
    for(len_t ip=0;ip<nShard;ip++){
        //if((rpPrevious[ip]>cutoffFrac*rp_initial[ip] && rp[ip]>cutoffFrac*rp_initial[ip]) && irp[ip]<nr){
        if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr){
            //rpdot[ip]=-lambda*pow(Tcold[irp[ip]]/2000.0,5.0/3.0)*pow(rp[ip]/0.002,4.0/3.0)*cbrt(ncold[irp[ip]]/1e20)/(4.0*M_PI*rp[ip]*rp[ip]*pelletDensity);
            //rpdot[ip]=-lambda*pow(Tcold[irp[ip]]/2000.0,5.0/3.0)*pow(rp[ip]/0.002,4.0/3.0)*cbrt(ncold[irp[ip]]/1e20)/(4.0*M_PI*(rp[ip]+cutoffFrac*rp_initial[ip])*(rp[ip]+cutoffFrac*rp_initial[ip])*pelletDensity);
            rpdot[ip]=-5/3*lambda*pow(Tcold[irp[ip]]/2000.0,5.0/3.0)*pow(1.0/0.002,4.0/3.0)*cbrt(ncold[irp[ip]]/1e20)/(4.0*M_PI*pelletDensity);
        }else
            rpdot[ip]=0;
    }
}

void SPIHandler::CalculateDepositionRate(){
    for(len_t ir=0;ir<nr;ir++){
        depositionRate[ir]=0;
        for(len_t ip=0;ip<nShard;ip++){
            //if((rpPrevious[ip]>cutoffFrac*rp_initial[ip] && rp[ip]>cutoffFrac*rp_initial[ip]) && irp[ip]<nr)
            if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr){
                // depositionRate[ir]+=-4.0*M_PI*rp[ip]*rp[ip]*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip];
                depositionRate[ir]+=-3.0/5.0*4.0*M_PI*pow(rp[ip],4.0/5.0)*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip];
            }
        }
    }
}

void SPIHandler::CalculateAdiabaticHeatAbsorbtionRateNGS(){
    for(len_t ir=0;ir<nr;ir++){
        heatAbsorbtionRate[ir]=0;
        for(len_t ip=0;ip<nShard;ip++){
            //if((rpPrevious[ip]>cutoffFrac*rp_initial[ip] && rp[ip]>cutoffFrac*rp_initial[ip]) && irp[ip]<nr)
            if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr)
                heatAbsorbtionRate[ir]+=-M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
        }
    }
}

void SPIHandler::CalculateTimeAveragedDeltaSourceLocal(real_t *timeAveragedDeltaSource){
    for(len_t ip=0;ip<nShard;ip++){
        rSourceMax=max(rCoordPPrevious[ip],rCoordPNext[ip]); 
        rSourceMin=min(rCoordPPrevious[ip],rCoordPNext[ip]);
        for(len_t ir=0;ir<nr;ir++){
            if((rGrid->GetR_f(ir)>rSourceMax || rGrid->GetR_f(ir+1)<rSourceMin)){
                timeAveragedDeltaSource[ir*nShard+ip]=0.0;
            }else{
                timeAveragedDeltaSource[ir*nShard+ip]=1.0/(rGrid->GetVpVol(ir)*VpVolNormFactor*(rSourceMax-rSourceMin))*
                                                     (min(rGrid->GetR_f(ir+1),rSourceMax)-max(rGrid->GetR_f(ir),rSourceMin))/rGrid->GetDr(ir);
            }
        }
    }
}

void SPIHandler::CalculateGaussianSourceLocal(real_t *gaussianSource){
    for(len_t ip=0;ip<nShard;ip++){
        for(len_t ir=0;ir<nr;ir++){
            gaussianSource[ir*nShard+ip]=((erf((rGrid->GetR_f(ir+1)-rCoordPNext[ip])/rCld[ip])-erf((rGrid->GetR_f(ir)-rCoordPNext[ip])/rCld[ip]))/2+
                                          (erf((-rGrid->GetR_f(ir+1)-rCoordPNext[ip])/rCld[ip])-erf((-rGrid->GetR_f(ir)-rCoordPNext[ip])/rCld[ip]))/2)/  //contribution on the "other" side of the magnetic axis
                                         (2*M_PI*M_PI*VpVolNormFactor*(rGrid->GetR_f(ir+1)*rGrid->GetR_f(ir+1)-rGrid->GetR_f(ir)*rGrid->GetR_f(ir)));
        }
    }
}

// This functionality should be moved to the radial grid generator-classes, and be optimised for every specific grid
void SPIHandler::CalculateIrp(){
    for(len_t ip=0;ip<nShard;ip++){
        for(len_t ir=0; ir<nr;ir++){
            if(rCoordPNext[ip]<rGrid->GetR_f(ir+1) && rCoordPNext[ip]>rGrid->GetR_f(ir)){
                irp[ip]=ir;
                break;
            }
            irp[ip]=nr;
        }
    }
}

void SPIHandler::CalculateRCld(){
    for(len_t ip=0;ip<nShard;ip++){
        if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT){
            rCld[ip]=rclPrescribedConstant;
        }else if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT){
            //rCld[ip]=10*rp[ip];// Very approximate, could be improved based on the Parks 2005 paper
            rCld[ip]=10*pow(rp[ip],3.0/5.0);// Very approximate, could be improved based on the Parks 2005 paper
        }
    }
}

real_t SPIHandler::CalculateLambda(real_t X){
    return (27.0837+tan(1.48709*X))/1000.0;
}

void SPIHandler::evaluatePartialContributionRpDot(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        evaluatePartialContributionRpDotNGS(jac, derivId, scaleFactor);
    }
}

void SPIHandler::evaluatePartialContributionDepositionRate(FVM::Matrix *jac, len_t derivId, real_t scaleFactor, real_t SPIMolarFraction, len_t rOffset){
    if((spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL || 
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE)||
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        evaluatePartialContributionDepositionRateNGS(jac, derivId, scaleFactor, SPIMolarFraction, rOffset);
    }
}

void SPIHandler::evaluatePartialContributionAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS ||
        spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN){
        evaluatePartialContributionAdiabaticHeatAbsorbtionRateNGS(jac, derivId, scaleFactor);
    }
}

void SPIHandler::evaluatePartialContributionRpDotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    if(derivId==id_rp){
        /*for(len_t ip=0;ip<nShard;ip++){
            //if((rpPrevious[ip]>cutoffFrac*rp_initial[ip] && rp[ip]>cutoffFrac*rp_initial[ip]) && irp[ip]<nr)
            if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr)
                jac->SetElement(ip,ip,scaleFactor* (-2.0/3.0)*rpdot[ip]/rp[ip]);
        }*/
    }else if(derivId==id_Tcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*5.0/3.0*rpdot[ip]/Tcold[irp[ip]]);
        }
    }else if(derivId==id_ncold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*rpdot[ip]/ncold[irp[ip]]);
        }
    }
}

void SPIHandler::evaluatePartialContributionDepositionRateNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor, real_t SPIMolarFraction, len_t rOffset){
    if(derivId==id_rp){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                //jac->SetElement(ir+rOffset,ip,-scaleFactor*SPIMolarFraction*4.0/3.0*4.0*M_PI*rp[ip]*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]);
                if(rp[ip]>cutoffFrac*rp_initial[ip])
                    jac->SetElement(ir+rOffset,ip,-scaleFactor*SPIMolarFraction*12.0/25.0*4.0*M_PI/pow(rp[ip],1.0/5.0)*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]);
            }
        }
    }else if(derivId==id_Tcold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr){
                    //jac->SetElement(ir+rOffset,irp[ip],-scaleFactor*SPIMolarFraction*5/3/Tcold[irp[ip]]*4*M_PI*rp[ip]*rp[ip]*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]);
                    jac->SetElement(ir+rOffset,irp[ip],-scaleFactor*SPIMolarFraction/Tcold[irp[ip]]*4.0*M_PI*pow(rp[ip],4.0/5.0)*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]);
                }
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr){
                    //jac->SetElement(ir+rOffset,irp[ip],-scaleFactor*SPIMolarFraction*1.0/3.0/ncold[irp[ip]]*4.0*M_PI*rp[ip]*rp[ip]*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]);
                    jac->SetElement(ir+rOffset,irp[ip],-scaleFactor*SPIMolarFraction*1.0/5.0/ncold[irp[ip]]*4.0*M_PI*pow(rp[ip],4.0/5.0)*rpdot[ip]*pelletDensity*Constants::N_Avogadro/pelletMolarMass*depositionProfilesAllShards[ir*nShard+ip]); 
               } 
           }
        }
    }
}

void SPIHandler::evaluatePartialContributionAdiabaticHeatAbsorbtionRateNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    if(derivId==id_rp){
        if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT){
            for(len_t ir=0;ir<nr;ir++){
                for(len_t ip=0;ip<nShard;ip++){
                    //if((rpPrevious[ip]>cutoffFrac*rp_initial[ip] && rp[ip]>cutoffFrac*rp_initial[ip]) && irp[ip]<nr)
                    if(rp[ip]>cutoffFrac*rp_initial[ip] && irp[ip]<nr){
                        //jac->SetElement(ir,ip,-scaleFactor*2.0/rp[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip]);
                        jac->SetElement(ir,ip,-scaleFactor*6.0/5.0/rp[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip]);
                    }
                }
            }
        }
    }else if(derivId==id_Tcold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(irp[ip]<nr)
                    jac->SetElement(ir,irp[ip],-scaleFactor*3.0/2.0*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*heatAbsorbtionProfilesAllShards[ir*nShard+ip]);
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(irp[ip]<nr)
                    jac->SetElement(ir,irp[ip],-scaleFactor*M_PI*rCld[ip]*rCld[ip]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip]);
            }
        }
    }
}
