/**
 * Implementation of a class that calculates and stores quantities related to the SPI shards
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include "DREAM/Equations/SPIHandler.hpp"

using namespace DREAM;
using namespace std;

/** 
 * List of molar masses and solid densities.
 * Needed to calculate the density of a mixed species pellet with molar fractions given.
 * This is in turn needed in the NGS formula since the ablation rate (from Parks TSDW 2017) is given in g/s
 * Isotope 0 means naturally occuring mix
 */
const len_t SPIHandler::nMolarMassList=3;
const len_t SPIHandler::ZMolarMassList[nMolarMassList]={1,1,10};
const len_t SPIHandler::isotopesMolarMassList[nMolarMassList]={2,0,0};// 0 means naturally occuring mix
const real_t SPIHandler::molarMassList[nMolarMassList]={0.0020141,0.001008,0.020183};// kg/mol

const len_t SPIHandler::nSolidDensityList=3;
const len_t SPIHandler::ZSolidDensityList[nSolidDensityList]={1,1,10};
const len_t SPIHandler::isotopesSolidDensityList[nSolidDensityList]={2,0,0};
const real_t SPIHandler::solidDensityList[nSolidDensityList]={205.9,86,1444};// kg/m^3

// Normalisation constants used in the NGS formula
const real_t T0=2000.0;// eV
const real_t n0=1e20;// m^{-3}
const real_t r0=0.002;// m

/**
 * Constructor
 */
SPIHandler::SPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t *Z, len_t *isotopes, const real_t *molarFraction, len_t NZ, 
    OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode,
    OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode,
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode,
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode,
    OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode, real_t VpVolNormFactor=1, real_t rclPrescribedConstant=0.01){

    // Get pointers to relevant objects
    this->rGrid=g->GetRadialGrid();
    this->unknowns=u;
    this->VpVolNormFactor=VpVolNormFactor;

    // Store settings
    this->spi_velocity_mode=spi_velocity_mode;
    this->spi_ablation_mode=spi_ablation_mode;
    this->spi_deposition_mode=spi_deposition_mode;
    this->spi_heat_absorbtion_mode=spi_heat_absorbtion_mode;
    this->spi_cloud_radius_mode=spi_cloud_radius_mode;

    // Set prescribed cloud radius (if any)
    if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT){
        this->rclPrescribedConstant=rclPrescribedConstant;
    }else{
        this->rclPrescribedConstant=0.0;
    }

    // Get unknown ID's
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    id_Yp = unknowns->GetUnknownID(OptionConstants::UQTY_Y_P);
    id_xp = unknowns->GetUnknownID(OptionConstants::UQTY_X_P);
    id_vp = unknowns->GetUnknownID(OptionConstants::UQTY_V_P);
    id_Wcold = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
  	    id_Whot = unknowns->GetUnknownID(OptionConstants::UQTY_W_HOT);
        id_qhot = unknowns->GetUnknownID(OptionConstants::UQTY_Q_HOT);
        id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    }
    // Get number of grid points and number of shards
    this->nr=rGrid->GetNr();
    this->nShard=unknowns->GetUnknown(id_Yp)->NumberOfMultiples();

    // Memory allocation
    AllocateQuantities();
    
    // Initialize rCoordPrevious to the radial coordinate at the plasma edge 
    // to use as a starting guess if rCoord must be solved for numerically
    for(len_t ip=0;ip<nShard;ip++)
        rCoordPPrevious[ip]=rGrid->GetR_f(nr-1);

    // Calculate pellet molar mass, molar volume and density
    real_t molarMass=0;
    real_t solidDensity=0;
    real_t *pelletDeuteriumFraction=new real_t[nShard];
    for(len_t ip=0;ip<nShard;ip++){
        pelletMolarMass[ip]=0;
        pelletMolarVolume[ip]=0;
    }
    
    len_t offset=0;
    for(len_t iZ=0;iZ<NZ;iZ++){
        if(molarFraction[offset]>=0){
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
            for(len_t ip=0;ip<nShard;ip++){
		        pelletMolarMass[ip]+=molarMass*molarFraction[offset+ip];
		        pelletMolarVolume[ip]+=molarMass/solidDensity*molarFraction[offset+ip];
		 
		        if(Z[iZ]==1 && isotopes[iZ]==2){
		            pelletDeuteriumFraction[ip]+=molarFraction[offset+ip];
		        }
            }
            offset+=nShard;
        }else {
        	offset+=1;
        }
    }
    
    for(len_t ip=0;ip<nShard;ip++){
		pelletDensity[ip]=pelletMolarMass[ip]/pelletMolarVolume[ip];
		
		/**
		* Evaluate the lambda factor that differs for different pellet compositions
		* It seems that the lambda implemented here is only valid for composite neon-deuterium pellets
		* but since the only reference for it is Parks 2017 TSDW presentation it is rather unclear.
		* Also note that lambda in Parks TSDW presentation is defined in terms of the molar fraction of D_2, 
		* while the input gives the molar fraction of D, hence the seemingly weird input argument.
		*/
		lambda[ip]=CalculateLambda(pelletDeuteriumFraction[ip]/2.0/(1.0-pelletDeuteriumFraction[ip]/2.0));
		
		if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS)
			NGSConstantFactor[ip]=5.0/3.0*lambda[ip]*pow(1.0/T0,5.0/3.0)*pow(1.0/r0,4.0/3.0)*cbrt(1.0/n0)/(4.0*M_PI*pelletDensity[ip]);
		else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS)
			NGSConstantFactor[ip]=5.0/3.0*pow(M_PI*Constants::me/256.0,1.0/6.0)*lambda[ip]*pow(1.0/(Constants::ec*T0),5.0/3.0)*pow(1.0/r0,4.0/3.0)*cbrt(1.0/n0)/(4.0*M_PI*pelletDensity[ip]);
	}
	delete [] pelletDeuteriumFraction;
}

/**
 * Destructor
 */
SPIHandler::~SPIHandler(){
   DeallocateQuantities();
}

/**
 * Allocate memory for arrays stored in this object
 */
void SPIHandler::AllocateQuantities(){
    DeallocateQuantities();

    Ypdot = new real_t[nShard];
    rCld = new real_t[nShard];
    depositionRate = new real_t[nr];
    depositionProfilesAllShards = new real_t[nr*nShard];
    heatAbsorbtionRate = new real_t[nr];
    heatAbsorbtionProfilesAllShards = new real_t[nr*nShard];
    rCoordPPrevious = new real_t[nShard];
    thetaCoordPPrevious = new real_t[nShard];
    rCoordPNext = new real_t[nShard];
    thetaCoordPNext = new real_t[nShard];
    irp = new len_t[nShard];
    gradRCartesian = new real_t[3];
    gradRCartesianPrevious = new real_t[3];
    qtot = new real_t[nr];
    Eeff = new real_t[nr];
    pelletMolarMass = new real_t[nShard];
    pelletMolarVolume = new real_t[nShard];
    pelletDensity = new real_t[nShard];
    lambda = new real_t[nShard];
    NGSConstantFactor = new real_t[nShard];
}

/**
 * Deallocate memory for arrays stored in this object
 */
void SPIHandler::DeallocateQuantities(){
    delete [] Ypdot;
    delete [] rCld;
    delete [] depositionRate;
    delete [] depositionProfilesAllShards;
    delete [] heatAbsorbtionRate;
    delete [] heatAbsorbtionProfilesAllShards;
    delete [] rCoordPPrevious;
    delete [] thetaCoordPPrevious;
    delete [] rCoordPNext;
    delete [] thetaCoordPNext;
    delete [] irp;
    delete [] gradRCartesian;
    delete [] gradRCartesianPrevious;
    delete [] qtot;
    delete [] Eeff;
    delete [] pelletMolarMass;
    delete [] pelletMolarVolume;
    delete [] pelletDensity;
    delete [] lambda;
    delete [] NGSConstantFactor;
}

/**
 * Rebuild this object
 * dt: current time step duration
 */
void SPIHandler::Rebuild(real_t dt){

    // Collect current data, and for some variables data from the previous time step
    xp=unknowns->GetUnknownData(id_xp);

    // Needed for the calculation of the time averaged delta function
    xpPrevious=unknowns->GetUnknownDataPrevious(id_xp);

    vp=unknowns->GetUnknownData(id_vp);
    ncold=unknowns->GetUnknownData(id_ncold);
    Tcold=unknowns->GetUnknownData(id_Tcold);
    Yp=unknowns->GetUnknownData(id_Yp);
    Wcold=unknowns->GetUnknownData(id_Wcold);
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
        Whot=unknowns->GetUnknownData(id_Whot);
        qhot=unknowns->GetUnknownData(id_qhot);
        ntot=unknowns->GetUnknownData(id_ntot);
    }

    // We use YpPrevious>0 as condition to keep the pellet terms active, 
    // to avoid making the functions discontinuous within a single time step
    YpPrevious=unknowns->GetUnknownDataPrevious(id_Yp);

    // We need the time step to calculate the transient factor in the deposition rate
    this->dt=dt;

    // Calculate current and previus radial coordinate from the cartesian coordinates used for the shard positions
    // (unlesa the shards do not have a velocity)
    if(spi_velocity_mode==OptionConstants::EQTERM_SPI_VELOCITY_MODE_PRESCRIBED){
        
        // We calculate the distance the shard travels during one time step, 
        // to be used as a length scale used to determine the tolerance 
        // when doing a numerical coordinate transformation
        real_t distP;
        
        // Note that at the first iteration, xp in the current time step will be equal to xpPrevious, unless xp is prescribed!
        for(len_t ip=0;ip<nShard;ip++){
            distP=sqrt((xp[3*ip]-xpPrevious[3*ip])*(xp[3*ip]-xpPrevious[3*ip])+
                       (xp[3*ip+1]-xpPrevious[3*ip+1])*(xp[3*ip+1]-xpPrevious[3*ip+1])+
                       (xp[3*ip+2]-xpPrevious[3*ip+2])*(xp[3*ip+2]-xpPrevious[3*ip+2]));
            rGrid->GetRThetaFromCartesian(&rCoordPPrevious[ip], &thetaCoordPPrevious[ip], xpPrevious[3*ip], xpPrevious[3*ip+1], xpPrevious[3*ip+2], distP, rCoordPPrevious[ip]);
            rGrid->GetRThetaFromCartesian(&rCoordPNext[ip], &thetaCoordPNext[ip], xp[3*ip], xp[3*ip+1], xp[3*ip+2], distP, rCoordPPrevious[ip]);
        }
    }else if(spi_velocity_mode==OptionConstants::EQTERM_SPI_VELOCITY_MODE_NONE){
        for(len_t ip=0;ip<nShard;ip++){
        	// If the shards do not move, we can not use the distance 
        	// the shards travel in one time step as a length scale to set the tolerance.
        	// Here we use a hardcoded length scale of 1 cm
            rGrid->GetRThetaFromCartesian(&rCoordPPrevious[ip], &thetaCoordPPrevious[ip], xpPrevious[3*ip], xpPrevious[3*ip+1], xpPrevious[3*ip+2], 0.01, rCoordPPrevious[ip]);
            rGrid->GetRThetaFromCartesian(&rCoordPNext[ip], &thetaCoordPNext[ip], xp[3*ip], xp[3*ip+1], xp[3*ip+2], 0.01, rCoordPPrevious[ip]);
        }
    }else {throw DREAMException("SPIHandler: unrecognized SPI shard velocity mode");}
    
    // Calculate the radial index of each shard
    CalculateIrp();
    
    // Calculate ablation rate (if any)
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        CalculateYpdotNGSParksTSDW();
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){ 
    
        for(len_t ir=0;ir<nr;ir++){
            // Total electron heat flux. 
            // The factor 1/4 is an approximate way to convert from the flux in all directions to the flux in only one direction
            qtot[ir]=(qhot[ir] + 4.0*sqrt(2.0/(M_PI*Constants::me))*ncold[ir]*pow(Constants::ec*Tcold[ir],3.0/2.0))/4.0;
            
            // Effective energy of incomming electrons
            Eeff[ir]=4.0/3.0*(Wcold[ir]+Whot[ir])/ntot[ir];
        }
        CalculateYpdotNGSParksTSDWKinetic();
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        for(len_t ip=0;ip<nShard;ip++)
            Ypdot[ip]=0;
    }else ifif(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGPS){
        throw NotImplementedException("SPIHandler: NGPS ablation is not yet supported");
    }else {throw DREAMException("SPIHandler: unrecognized SPI shard ablation mode");}

    // Calculate radius of the neutral cloud (if any)
    if(spi_cloud_radius_mode!=OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_NEGLECT)
        CalculateRCld();

    // Calculate deposition (if any)
    if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL){
        CalculateTimeAveragedDeltaSourceLocal(depositionProfilesAllShards);

    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){
        CalculateTimeAveragedDeltaSourceLocal(depositionProfilesAllShards);

        // Shift the deposition profile to the last grid cell (before the current one) to avoid "self-dilution"
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


    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        CalculateGaussianSourceLocal(depositionProfilesAllShards);

    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT){
        for(len_t ir=0;ir<nr;ir++)
            depositionRate[ir]=0;
    }else {throw DREAMException("SPIHandler: unrecognized SPI material deposition mode");}

    // Calculate heat absorbtion
    if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS){
        CalculateTimeAveragedDeltaSourceLocal(heatAbsorbtionProfilesAllShards);
        CalculateAdiabaticHeatAbsorbtionRateMaxwellian();

    }else if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN){
        CalculateGaussianSourceLocal(heatAbsorbtionProfilesAllShards);
        CalculateAdiabaticHeatAbsorbtionRateMaxwellian();

    }else if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT){
        for(len_t ir=0;ir<nr;ir++)
            heatAbsorbtionRate[ir]=0;
    }else {throw DREAMException("SPIHandler: unrecognized SPI heat absorbtion mode");}
}

/**
 * Deposition rate according to NGS formula from Parks TSDW presentation 2017
 */
void SPIHandler::CalculateYpdotNGSParksTSDW(){
    for(len_t ip=0;ip<nShard;ip++){
        if(YpPrevious[ip]>0 && irp[ip]<nr){
            Ypdot[ip]=-NGSConstantFactor[ip]*pow(Tcold[irp[ip]],5.0/3.0)*cbrt(ncold[irp[ip]]);
        }else
            Ypdot[ip]=0;
    }
}

/**
 * Deposition rate according to NGS formula from Parks TSDW presentation 2017,
 * including contribution from a kinetically treated species. The kinetically treated
 * species is taken into account by expressing the NGS formula in terms of the
 * heat flux and effective energy of the incomming electrons, and then calculating the
 * total value of these for the hot and cold population combined.
 */
void SPIHandler::CalculateYpdotNGSParksTSDWKinetic(){
    for(len_t ip=0;ip<nShard;ip++){
        if(YpPrevious[ip]>0 && irp[ip]<nr){
            Ypdot[ip]=-NGSConstantFactor[ip]*pow(qtot[irp[ip]],1.0/3.0)*pow(Eeff[irp[ip]],7.0/6.0);
        }else
            Ypdot[ip]=0;
    }
}

/**
 * Calculate deposition corresponding to the ablation, with a density conserving discretisation
 */
real_t *SPIHandler::CalculateDepositionRate(real_t *SPIMolarFraction){
    for(len_t ir=0;ir<nr;ir++){
        depositionRate[ir]=0;
        for(len_t ip=0;ip<nShard;ip++){
            if(YpPrevious[ip]>0 && irp[ip]<nr){
                depositionRate[ir]+=-SPIMolarFraction[ip]*4.0*M_PI*(Yp[ip]/abs(Yp[ip])*pow(abs(Yp[ip]),9.0/5.0)-pow(YpPrevious[ip],9.0/5.0))/3.0/pelletMolarVolume[ip]*Constants::N_Avogadro/dt*depositionProfilesAllShards[ir*nShard+ip];
            }
        }
    }
    return depositionRate;
}

/**
 * Calculate the total heat flux going into the pellet cloud assuming Maxwellian distribution for the insident electrons
 */
void SPIHandler::CalculateAdiabaticHeatAbsorbtionRateMaxwellian(){
    for(len_t ir=0;ir<nr;ir++){
        heatAbsorbtionRate[ir]=0;
        for(len_t ip=0;ip<nShard;ip++){
            if(YpPrevious[ip]>0 && irp[ip]<nr)
                heatAbsorbtionRate[ir]+=-M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
        }
    }
}

/**
 * Calculate a delta function averaged over the current time step (which gives a "box"-function) 
 * and grid cell volume (splits the "box" between the cells passed during the time step
 * Derivation available on github in SPIDeltaSoure.pdf
 */
void SPIHandler::CalculateTimeAveragedDeltaSourceLocal(real_t *timeAveragedDeltaSource){
    for(len_t ip=0;ip<nShard;ip++){

        // Reset the deposition profile
        for(len_t ir=0;ir<nr;ir++)
            timeAveragedDeltaSource[ir*nShard+ip]=0.0;

        // Find out if the shard has passed its turning point (where the radial coordinate goes from decreasing to increasing)
        // If so, we split the averaging in two steps, before and after the turning point (where the analythical expression used for the averaging breaks down)
        // We determine wether the turning point has been passed by checking if there is a sign change of the 
        // projection of the direction vector for the shard motion on the gradient of the radial coordinate before and after the time step
        nSplit=1;
        iSplit=0;
        turningPointPassed=false;
        rGrid->GetGradRCartesian(gradRCartesian,rCoordPNext[ip],thetaCoordPNext[ip]);
        rGrid->GetGradRCartesian(gradRCartesianPrevious,rCoordPPrevious[ip],thetaCoordPPrevious[ip]);
        turningPointPassed=((gradRCartesian[0]*(xp[3*ip]-xpPrevious[3*ip])+
            gradRCartesian[1]*(xp[3*ip+1]-xpPrevious[3*ip+1])+
            gradRCartesian[2]*(xp[3*ip+2]-xpPrevious[3*ip+2])) *
            (gradRCartesianPrevious[0]*(xp[3*ip]-xpPrevious[3*ip])+
            gradRCartesianPrevious[1]*(xp[3*ip+1]-xpPrevious[3*ip+1])+
            gradRCartesianPrevious[2]*(xp[3*ip+2]-xpPrevious[3*ip+2])))<0;

        while(iSplit<nSplit){
            if(turningPointPassed){
                if(iSplit==1){
                    nSplit=2;
                    rCoordPClosestApproach=rGrid->FindClosestApproach(xp[3*ip],xp[3*ip+1],xp[3*ip+2],xpPrevious[3*ip],xpPrevious[3*ip+1],xpPrevious[3*ip+2]);
                    rSourceMax=max(rCoordPPrevious[ip],rCoordPClosestApproach); 
                    rSourceMin=min(rCoordPPrevious[ip],rCoordPClosestApproach);
                }else{
                    rSourceMax=max(rCoordPClosestApproach,rCoordPNext[ip]); 
                    rSourceMin=min(rCoordPClosestApproach,rCoordPNext[ip]);
                }
            }else{
                rSourceMax=max(rCoordPPrevious[ip],rCoordPNext[ip]); 
                rSourceMin=min(rCoordPPrevious[ip],rCoordPNext[ip]);
            }
    
            for(len_t ir=0;ir<nr;ir++){
                if(!(rGrid->GetR_f(ir)>rSourceMax || rGrid->GetR_f(ir+1)<rSourceMin)){
                    timeAveragedDeltaSource[ir*nShard+ip]+=1.0/(rGrid->GetVpVol(ir)*VpVolNormFactor*(rSourceMax-rSourceMin))*
                                                         (min(rGrid->GetR_f(ir+1),rSourceMax)-max(rGrid->GetR_f(ir),rSourceMin))/rGrid->GetDr(ir);
                }
            }
            iSplit++;
        }
    }
}

/**
 * Calculates a gaussian deposition profile with 1/e length scale equal to the shards' cloud radii
 * NOTE: not time averaged, so be careful with using time steps long enough to allow the shards 
 * to travel lengths comparable to the cloud radius! Also, this profile is gaussian in the radial coordinate,
 * not a 2D-gaussian in the poloidal plane!
 */
void SPIHandler::CalculateGaussianSourceLocal(real_t *gaussianSource){
    for(len_t ip=0;ip<nShard;ip++){
        for(len_t ir=0;ir<nr;ir++){
            gaussianSource[ir*nShard+ip]=((erf((rGrid->GetR_f(ir+1)-rCoordPNext[ip])/rCld[ip])-erf((rGrid->GetR_f(ir)-rCoordPNext[ip])/rCld[ip]))/2+
                                          (erf((-rGrid->GetR_f(ir+1)-rCoordPNext[ip])/rCld[ip])-erf((-rGrid->GetR_f(ir)-rCoordPNext[ip])/rCld[ip]))/2)/  //contribution on the "other" side of the magnetic axis
                                         (2*M_PI*M_PI*VpVolNormFactor*(rGrid->GetR_f(ir+1)*rGrid->GetR_f(ir+1)-rGrid->GetR_f(ir)*rGrid->GetR_f(ir)));
        }
    }
}

/**
 * General function to find the grid cell indexes corresponding to the shard positions
 * This functionality could be moved to the radial grid generator-classes, and be optimised for every specific grid
 */
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

/**
 * Calculates the shards' cloud radius (no good way to do this self-consistently yet!)
 */
void SPIHandler::CalculateRCld(){
    for(len_t ip=0;ip<nShard;ip++){
        if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT){
            rCld[ip]=rclPrescribedConstant;
        }else if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT){

            // Very approximate. 
            // Could be improved based on the Parks 2005 paper, 
            // but the scaling given there does not agree with more advanced studies (eg Legnyel et al, NF 1999)
            rCld[ip]=10*pow(Yp[ip],3.0/5.0);
        }
    }
}

/**
 * Calculate lambda factor that differs between different pellet compositions
 * according to Parks 2017 TSDW presentation
 *
 * X: fraction of D2
 */
real_t SPIHandler::CalculateLambda(real_t X){
    return (27.0837+tan(1.48709*X))/1000.0;
}

/**
 * Wrapper function for calculation of partial derivatives of the ablation rate, depending on the settings
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */ 
void SPIHandler::evaluatePartialContributionYpdot(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        evaluatePartialContributionYpdotNGS(jac, derivId, scaleFactor);
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
        evaluatePartialContributionYpdotNGSKinetic(jac, derivId, scaleFactor);
    }
}

/**
 * Wrapper function for calculation of partial derivatives of the deposition rate, depending on the settings
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS and to apply weights (eg ionization energies and/or molar fractions)
 */
void SPIHandler::evaluatePartialContributionDepositionRate(FVM::Matrix *jac, len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset){
    if((spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL || 
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE)||
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        evaluatePartialContributionDepositionRateDensCons(jac, derivId, scaleFactor, SPIMolarFraction, rOffset);
    }
}

/**
 * Wrapper function for calculation of partial derivatives of the heat absorbtion rate, depending on the settings
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
void SPIHandler::evaluatePartialContributionAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS ||
        spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN){
        evaluatePartialContributionAdiabaticHeatAbsorbtionRateMaxwellian(jac, derivId, scaleFactor);
    }
}

/**
 * Sets jacobian block contribution from the ablation rate,
 * for the NGS ablation rate from Parks 2017 TSDW presentation
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
void SPIHandler::evaluatePartialContributionYpdotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    if(derivId==id_Tcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*5.0/3.0*Ypdot[ip]/Tcold[irp[ip]]);
        }
    }else if(derivId==id_ncold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/ncold[irp[ip]]);
        }
    }
}

/**
 * Sets jacobian block contribution from the ablation rate,
 * for the NGS ablation rate from Parks 2017 TSDW presentation,
 * with a kinetic hot population included
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
void SPIHandler::evaluatePartialContributionYpdotNGSKinetic(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    if(derivId==id_Tcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)        		
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]*3.0/2.0*(qtot[irp[ip]]-qhot[irp[ip]])/Tcold[irp[ip]]);
        }
    }else if(derivId==id_ncold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]*(qtot[irp[ip]]-qhot[irp[ip]])/ncold[irp[ip]]);
        }
    }else if(derivId==id_Whot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*7.0/6.0*Ypdot[ip]/(Whot[irp[ip]]+Wcold[irp[ip]]));
        }
    }else if(derivId==id_qhot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]);
        }
    }else if(derivId==id_ntot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], -scaleFactor*7.0/6.0*Ypdot[ip]/ntot[irp[ip]]);
        }
    }else if(derivId==id_Wcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr)
                jac->SetElement(ip,irp[ip], scaleFactor*7.0/6.0*Ypdot[ip]/(Whot[irp[ip]]+Wcold[irp[ip]]));
        }
    }
}

/**
 * Sets jacobian block contribution from the deposition rate,
 * for the density conserving discretisation
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS and to apply weights (eg ionization energies and/or molar fractions)
 * SPIMolarFraction: molar fraction of the pellet consisting of the currently considered species
 * rOffset: offset for the currently considered species in the vector with ion densities at different radial indexes
 */
void SPIHandler::evaluatePartialContributionDepositionRateDensCons(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset){
    if(derivId==id_Yp){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(YpPrevious[ip]>0)
                    jac->SetElement(ir+rOffset,ip,-scaleFactor[ir]*SPIMolarFraction[ip]*12.0/5.0*M_PI*pow(abs(Yp[ip]),4.0/5.0)/pelletMolarVolume[ip]*Constants::N_Avogadro/dt*depositionProfilesAllShards[ir*nShard+ip]);
            }
        }
    }
}

/**
 * Sets jacobian block contribution from the heat absorbtion rate,
 * assuming Maxwellian distribution for the insident electrons
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
void SPIHandler::evaluatePartialContributionAdiabaticHeatAbsorbtionRateMaxwellian(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    if(derivId==id_Yp){
        if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT){
            for(len_t ir=0;ir<nr;ir++){
                for(len_t ip=0;ip<nShard;ip++){
                    if(YpPrevious[ip]>0 && irp[ip]<nr){
                        jac->SetElement(ir,ip,-scaleFactor*6.0/5.0/Yp[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]*heatAbsorbtionProfilesAllShards[ir*nShard+ip]);
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
