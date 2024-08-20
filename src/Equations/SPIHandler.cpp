/**
 * Implementation of a class that calculates and stores quantities related to the SPI shards
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include "DREAM/Equations/SPIHandler.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include "DREAM/IO.hpp"

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

const real_t qe = Constants::ec;//1.60217662e-19;// C
const real_t me = Constants::me;//9.10938356e-31;// kg
const real_t gamma_e = 1;//Adiabatic constant of electrons
const real_t gamma_i = 3;//Adiabatic constant of ions
const real_t Zavg0Drift = 1;// Assume all ion species are singly ionized directly after the neutral phase and while the cloud detaches from the pellet
const real_t N_Avogadro = Constants::N_Avogadro;//6.02214076e23

/**
 * Constructor
 */
SPIHandler::SPIHandler(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t *Z, len_t *isotopes, const real_t *molarFraction, len_t NZ, 
    OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode,
    OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode,
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode,
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode,
    OptionConstants::eqterm_spi_cloud_radius_mode spi_cloud_radius_mode,
    OptionConstants::eqterm_spi_magnetic_field_dependence_mode spi_magnetic_field_dependence_mode, 
    OptionConstants::eqterm_spi_shift_mode spi_shift_mode, 
    real_t *TDrift, real_t T0Drift, real_t DeltaYDrift, real_t RmDrift, real_t *ZavgDriftArray,
	len_t nZavgDrift, len_t *ZsDrift, len_t *isotopesDrift,
    real_t VpVolNormFactor=1, real_t rclPrescribedConstant=0.01, const int_t *nbrShiftGridCell=nullptr){

    // Get pointers to relevant objects
    this->rGrid=g->GetRadialGrid();
    this->unknowns=u;
    this->VpVolNormFactor=VpVolNormFactor;

	// Get the major radius, to be used to properly normalize VpVol
	real_t R0 = this->rGrid->GetR0();
    rf = this->rf;
    qBgDrift = 1; //TODO Update if a self-consistent calculation of the q-profile becomes available

	// If R0 is infinite, i.e. toroidicity is not included in the simulation,
	// we can not use R0 from the radial grid of this simulation to calculate 
	// the size of the flux surfaces. The corresponding factor correcting the 
	// size of the flux surfaces must instead be included directly in the 
	// VpVolNormFactor. We also need the major radius for the drift displacement
	// calculation, so we also set this here if the R0 from the radial grid is finite
	// and no other value is specified by the user.
	if(!isinf(R0)){
	    this->VpVolNormFactor*=R0;
        if(RmDrift==-1)// RmDrift==-1 meeans that no value has been given by the user
            RmDrift=R0;// RmDrift is the major radius used for the drift calculation
    }else if(RmDrift==-1 && spi_shift_mode!=OptionConstants::EQTERM_SPI_SHIFT_MODE_NEGLECT)
        throw DREAMException("SPIHandler: The drift model requires a finite major radius.");

    // Store settings
    this->spi_velocity_mode=spi_velocity_mode;
    this->spi_ablation_mode=spi_ablation_mode;
    this->spi_deposition_mode=spi_deposition_mode;
    this->spi_heat_absorbtion_mode=spi_heat_absorbtion_mode;
    this->spi_cloud_radius_mode=spi_cloud_radius_mode;
    this->spi_magnetic_field_dependence_mode=spi_magnetic_field_dependence_mode;
    this->spi_shift_mode=spi_shift_mode;

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
    id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_Wcold = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
  	    id_Whot = unknowns->GetUnknownID(OptionConstants::UQTY_W_HOT);
        id_qhot = unknowns->GetUnknownID(OptionConstants::UQTY_Q_HOT);
        id_ntot = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    }
    // Get number of grid points and number of shards
    this->nr=rGrid->GetNr();
    this->nShard=unknowns->GetUnknown(id_Yp)->NumberOfMultiples();
    
    this->NZ=NZ;

    // Memory allocation
    AllocateQuantities();
    
    // Ablation cloud quantities
    for(len_t ip=0;ip<nShard;ip++)
        this->TDrift[ip]=TDrift[ip];
    this->T0Drift=T0Drift;
    this->DeltaYDrift=DeltaYDrift;
    this->RmDrift=RmDrift;
    this->ZavgDriftArray=ZavgDriftArray;
    
    
    // Initialize rCoordPrevious to the radial coordinate at the plasma edge 
    // to use as a starting guess if rCoord must be solved for numerically
    for(len_t ip=0;ip<nShard;ip++)
        rCoordPPrevious[ip]=rGrid->GetR_f(nr-1);

    // Calculate pellet molar mass, molar volume, density and average charge of the ablation cloud
    // The later is calculated based on a user-given lookup table (ZavgArray) for the average charge
    // of the relevant species (specified in Zs and isotopesDrift), as the ADAS rates are not valid 
    // for the conditions inside the ablation cloud
    real_t molarMass=0;
    real_t solidDensity=0;
    real_t ZavgDriftList=0;
    for(len_t ip=0;ip<nShard;ip++){
        pelletMolarMass[ip]=0;
        pelletMolarVolume[ip]=0;
        pelletDeuteriumFraction[ip]=0;
        ZavgDrift[ip]=0;
    }
    
    len_t offset=0;
    len_t counter1=0;
    len_t counter2=0;
    for(len_t iZ=0;iZ<NZ;iZ++){
        if(molarFraction[offset]>=0){
            for(len_t i=0;i<nMolarMassList;i++){
                if(Z[iZ]==ZMolarMassList[i] && isotopes[iZ]==isotopesMolarMassList[i]){
                    molarMass=molarMassList[i];
                }
                else{counter1++;}
            }
            if(counter1==nMolarMassList){throw DREAMException("SPIHandler: Pellet type is not recognized. Currently only neon and deuterium pellets are supported. To support other types fill in the material data in src/Equations/SPIHandler.cpp and py/DREAM/Settings/Equations/SPI.py.");}
            counter1=0;
            
            for(len_t i=0;i<nSolidDensityList;i++){
                if(Z[iZ]==ZSolidDensityList[i] && isotopes[iZ]==isotopesSolidDensityList[i]){
                    solidDensity=solidDensityList[i];
                }
            }
            for(len_t i=0;i<nZavgDrift;i++){
                if(Z[iZ]==ZsDrift[i] && isotopes[iZ]==isotopesDrift[i]){
                    ZavgDriftList=ZavgDriftArray[i];
                }
            }
            for(len_t ip=0;ip<nShard;ip++){
		        pelletMolarMass[ip]+=molarMass*molarFraction[offset+ip];
		        pelletMolarVolume[ip]+=molarMass/solidDensity*molarFraction[offset+ip];
                
		        if(Z[iZ]==1 && isotopes[iZ]==2)
		            pelletDeuteriumFraction[ip]+=molarFraction[offset+ip];

                ZavgDrift[ip]+=ZavgDriftList*molarFraction[offset+ip];
                if (!((Z[iZ]==1 && isotopes[iZ]==2) || Z[iZ]==10) && counter2==0){
                    DREAM::IO::PrintWarning(DREAM::IO::WARNING_ABLATION_RATE_NOT_VALID, "SPIHandler: The currently available ablation rate is only derived for deuterium and neon mixtures. Results may be inaccurate");
                    counter2++;
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
	
	// Set the presdcribed number of grid cells to shift the deposition for every shard (if any)
	if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE)
		for(len_t ip=0;ip<nShard;ip++)
	        this->nbrShiftGridCellPrescribed[ip] = 1;
	else if(spi_shift_mode==OptionConstants::EQTERM_SPI_SHIFT_MODE_PRESCRIBED) {
	    for(len_t ip=0;ip<nShard;ip++){
	        this->nbrShiftGridCellPrescribed[ip] = nbrShiftGridCell[ip];
        }
	
		delete [] nbrShiftGridCell;
	}

	delete [] ZsDrift;
	delete [] isotopesDrift;
}

/**
 * Destructor
 */
SPIHandler::~SPIHandler(){
	DeallocateQuantities();

	if (this->ZavgDriftArray != nullptr)
		delete [] this->ZavgDriftArray;
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
    phiCoordPPrevious = new real_t[nShard];
    rCoordPNext = new real_t[nShard];
    thetaCoordPNext = new real_t[nShard];
    phiCoordPNext = new real_t[nShard];
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
    nbrShiftGridCell = new int_t[nShard];
    nbrShiftGridCellPrescribed = new int_t[nShard];
    TDrift = new real_t[nShard];
    pelletDeuteriumFraction=new real_t[nShard];
    rp=new real_t[nShard];
    rpdot=new real_t[nShard];
    shift_r=new real_t[nShard];
    shift_store=new real_t[nShard];
    YpdotPrevious=new real_t[nShard];
    ZavgDrift=new real_t[nShard];
    plasmoidAbsorbtionFactor=new real_t[nShard];
    cosThetaDrift = new real_t[nShard];
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
    delete [] phiCoordPPrevious;
    delete [] rCoordPNext;
    delete [] thetaCoordPNext;
    delete [] phiCoordPNext;
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
    delete [] nbrShiftGridCell;
    delete [] nbrShiftGridCellPrescribed;
    delete [] TDrift;
    delete [] pelletDeuteriumFraction;
    delete [] rp;
    delete [] rpdot;
    delete [] shift_r;
    delete [] shift_store;
    delete [] YpdotPrevious;
    delete [] ZavgDrift;
    delete [] plasmoidAbsorbtionFactor;
    delete [] cosThetaDrift;
}
/**
* Calculates the radius and ablation of each shard
*/
void SPIHandler::YpConversion(len_t ip){
    rp[ip] = pow(YpPrevious[ip], 3.0/5.0);
    rpdot[ip] = 3.0/5.0 * pow(rp[ip], -2.0/3.0) * YpdotPrevious[ip];
}
/**
 * Below are functions needed to compute the drift displacement of an ionised cloud
 * succeeding the pellet injection. The analytical model is derived in doi:10.1017/S0022377823000466
 */

 // Stores data about the surroundings of a shard
void SPIHandler::AssignShardSpecificDriftParameters(len_t ip){
    v0Drift = -vp[3*ip];
    neBgDrift = ncoldPrevious[irp[ip]];
    TeBgDrift = TcoldPrevious[irp[ip]];
    BBgDrift = sqrt(this->rGrid->GetFSA_B2(irp[ip])) * rGrid->GetBmin(irp[ip]);
    sigmaBgDrift = rf->GetElectricConductivity(irp[ip]);
    niBgDrift = 0;
    for(len_t iZ=0;iZ<NZ;iZ++){
        niBgDrift += rf->GetIonHandler()->GetTotalIonDensity(irp[ip], iZ);
    }
}

/**
 * Computes relevant quantities for the drift displacement and heat absorbtion of a shard
 * Zavg    : Total average charge inside the cloud
 * CST     : Sound speed inside the cloud during the majority of the drift
 * CST0    : Sound speed inside the cloud directly after the neutral phase
 * G       : Mass ablation rate
 * n_0     : Density directly after the neutral phase
 * a0      : Initial cloud acceleration
 * t_detach: Time it takes to move away from the pellet (shard)
 * Lc      : Initial length
 * n       : Line integrated density
 * v_lab   : Initial radial drift velocity in the lab frame
 * Reff    : Effective resistance to ohmic currents exiting the cloud parallell to the field lines
 */
void SPIHandler::AssignDriftComputationParameters(len_t ip){
    CSTDrift = sqrt((gamma_e*ZavgDrift[ip] + gamma_i) * qe * TDrift[ip]/(pelletMolarMass[ip]/N_Avogadro));
    CST0Drift = sqrt((gamma_e*ZavgDrift[ip] + gamma_i) * qe * T0Drift/(pelletMolarMass[ip]/N_Avogadro));
    G = -4 * M_PI * pelletDensity[ip] * rp[ip] * rp[ip] * rpdot[ip];
    n0Drift = (1 + Zavg0Drift)*G/(2 * M_PI * DeltaYDrift * DeltaYDrift * pelletMolarMass[ip]/N_Avogadro * CST0Drift);
    a0Drift = ((1 + Zavg0Drift)*qe*T0Drift/(pelletMolarMass[ip]/N_Avogadro*RmDrift));
    tDetachDrift = -v0Drift/a0Drift + sqrt(v0Drift*v0Drift/(a0Drift*a0Drift) + 2*DeltaYDrift/a0Drift);
    LcInitDrift = 2 * CST0Drift*tDetachDrift;
    nBarDrift = n0Drift * LcInitDrift; 
    vLabInitDrift = a0Drift * tDetachDrift;
    ReffBgDrift = -2*M_PI*M_PI*RmDrift*this->rGrid->GetMinorRadius()/(sigmaBgDrift*DeltaYDrift*DeltaYDrift*DeltaYDrift*log(DeltaYDrift/(this->rGrid->GetMinorRadius()*M_PI)));
    
       
    // Calculate the fraction of the heat flux being absorbed in the cloud, using eq. 4.77 in Nicos MSc thesis
    real_t x = rf->GetElectronCollisionTimeThermal(irp[ip])*neBgDrift/n0Drift*(1+Zavg0Drift)*sqrt(2*TeBgDrift*qe/me)/LcInitDrift;// Ratio of mean free path and cloud length
    plasmoidAbsorbtionFactor[ip] = 1-1/(x*x)*(exp(-1/sqrt(x))*(-0.5*sqrt(x)+0.5*x+x*sqrt(x)+x*x)-0.5*expint(-1/sqrt(x)));

}

/**
 * Computes the characteristic quantities of time for the drift from a shard
 * tAccDrift : Acceleration time scale
 * tPolDrift : Poloidal rotation time scale
 * tPeDrift  : Pressure equilibration time scale
 * tExpDrift : Expansion time scale
 * tPolDriftPrime: Normalized tPolDrift
 * tPeDriftPrime : Normalized tPeDrift
 * tExpDriftPrime: Normalized tExpDrift
 */
void SPIHandler::AssignDriftTimeParameters(len_t ip){
    tAccDrift = nBarDrift/(1+ZavgDrift[ip])*pelletMolarMass[ip]/N_Avogadro*ReffBgDrift/(BBgDrift*BBgDrift);
    tPolDrift = qBgDrift*RmDrift/CSTDrift;
    tPeDrift = qe*TDrift[ip]*nBarDrift/(2*CSTDrift*(niBgDrift+neBgDrift)*qe*TeBgDrift);
    tExpDrift = LcInitDrift/(2*CSTDrift);
    tPolDriftPrime = tPolDrift/tAccDrift;
    tPeDriftPrime = tPeDrift/tAccDrift;
    tExpDriftPrime = tExpDrift/tAccDrift;
}


real_t SPIHandler::Integrand(real_t x, void *p){
    struct integrand_struct *params = (struct integrand_struct *)p;
    real_t temp = (params->a*cos(x) + x*sin(x))/((params->a)*(params->a)+x*x);
    return temp;
}
/**
 * Function to evaluate the difference of two complex exponential integrals in equation (A2) in doi:10.1017/S0022377823000466
 */
real_t SPIHandler::Epsiloni(real_t a, real_t b){
    real_t sum, error;
    if (a!=0){
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
        gsl_function F;
        integrand_struct paramstruct = {a};
        F.function = &Integrand;
        F.params = &paramstruct;
        if (b>0){
            gsl_integration_qags(&F, 0, b, 0, 1e-7, 1000, workspace, &sum, &error);
        }else{
            gsl_integration_qags(&F, b, 0, 0, 1e-7, 1000, workspace, &sum, &error);
            sum = -sum;
        }
        gsl_integration_workspace_free(workspace);
    }else if(b>=0){
        sum = gsl_sf_Si(b) + 0.5*M_PI;
    }else{
        sum = gsl_sf_Si(b) - 0.5*M_PI;
    }
    return sum;

}


real_t SPIHandler::BisFunction(real_t t_prime){
    return t_prime + tExpDriftPrime;
}

//Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveFirstRow(real_t t_prime){//- at cos in article
    real_t t_bis = BisFunction(t_prime);
    real_t term1 = Epsiloni(t_bis, t_bis/tPolDriftPrime);
    real_t numerator = (sin(t_bis/tPolDriftPrime) + 1/tPolDriftPrime * cos(t_bis/tPolDriftPrime));
    real_t denominator =  tPeDriftPrime * (1 + 1/(tPolDriftPrime * tPolDriftPrime));
    real_t term2 = numerator/denominator;
    return term1 - term2;
}

// Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveSecondRow(real_t t_prime){
    real_t t_bis = BisFunction(t_prime);
    real_t term1 = Epsiloni(0, t_bis/tPolDriftPrime);
    real_t term2 = Epsiloni(t_bis, t_bis/tPolDriftPrime);
    real_t result = term1 - term2;
    return result;
}

// Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveThirdRow(real_t t_prime){
    real_t t_bis = BisFunction(t_prime);
    real_t term1 = tPolDriftPrime * cos(t_bis/tPolDriftPrime);
    real_t term2 = sin(t_bis/tPolDriftPrime);
    return term1 + term2;
}

real_t SPIHandler::FirstRow(){
    real_t term1 = PrimitiveFirstRow(tPeDriftPrime);
    real_t term2 = PrimitiveFirstRow(0);
    return term1 - term2;
}

real_t SPIHandler::SecondRow(){
    real_t term1 = PrimitiveSecondRow(tPeDriftPrime);
    real_t term2 = PrimitiveSecondRow(0);
    return term1 - term2;
}

real_t SPIHandler::ThirdRow(){
    real_t factor = 1/tPeDriftPrime * 1/(1 + 1/(tPolDriftPrime*tPolDriftPrime));
    real_t term1 = PrimitiveThirdRow(tPeDriftPrime);
    real_t term2 = PrimitiveThirdRow(0);
    return factor * (term1 - term2);
}

// Function to collect all terms to evaluate equation A4
real_t SPIHandler::Deltar(len_t ip){
    real_t first = FirstRow();
    real_t second = SecondRow();
    real_t third = ThirdRow();
    real_t term1 = vLabInitDrift * tAccDrift;
    real_t factor = (1+ZavgDrift[ip])*2*qe*TDrift[ip]*qBgDrift/(CSTDrift*pelletMolarMass[ip]/N_Avogadro)*tAccDrift;

	real_t DeltaR = (term1 + factor * (first + second + third));
	real_t Rp = hypot(this->xp[ip*3+0], this->xp[ip*3+2]) + DeltaR;
	real_t yp = this->xp[ip*3+1];

	// Locate flux surface to which the pellet drifts
	real_t r=0, theta=0, phi=0;
	this->rGrid->GetRThetaPhiFromCartesian(&r, &theta, &phi, Rp, yp, 0, 0.01, Rp);

	real_t Deltar = r - this->rCoordPNext[ip];

    return Deltar;
}

/**
 * Rebuild this object
 * dt: current time step duration
 */
void SPIHandler::Rebuild(real_t dt, real_t t){

    // Collect current data, and for some variables data from the previous time step
    xp=unknowns->GetUnknownData(id_xp);

    // Needed for the calculation of the time averaged delta function
    xpPrevious=unknowns->GetUnknownDataPrevious(id_xp);

    vp=unknowns->GetUnknownData(id_vp);
    ncold=unknowns->GetUnknownData(id_ncold);
    Tcold=unknowns->GetUnknownData(id_Tcold);
    Yp=unknowns->GetUnknownData(id_Yp);
    Wcold=unknowns->GetUnknownData(id_Wcold);
    ni=unknowns->GetUnknownData(id_ni);
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
        Whot=unknowns->GetUnknownData(id_Whot);
        qhot=unknowns->GetUnknownData(id_qhot);
        ntot=unknowns->GetUnknownData(id_ntot);
    }

    // We use YpPrevious>0 as condition to keep the pellet terms active, 
    // to avoid making the functions discontinuous within a single time step
    YpPrevious=unknowns->GetUnknownDataPrevious(id_Yp);
    ncoldPrevious=unknowns->GetUnknownDataPrevious(id_ncold);
    TcoldPrevious=unknowns->GetUnknownDataPrevious(id_Tcold);

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
                       
            // If distP is zero (or at least very small), perhaps because some shards have not started to move yet
            // it is not suitable to use as a length scale to set the tolerance
            // Here we instead use a hardcoded length scale of 1 cm
            if(distP<1e-20)
            	distP=0.01;
            	
            if(t!=t_old)
                rGrid->GetRThetaPhiFromCartesian(&rCoordPPrevious[ip], &thetaCoordPPrevious[ip], &phiCoordPPrevious[ip], xpPrevious[3*ip], xpPrevious[3*ip+1], xpPrevious[3*ip+2], distP, rCoordPPrevious[ip]);
            rGrid->GetRThetaPhiFromCartesian(&rCoordPNext[ip], &thetaCoordPNext[ip], &phiCoordPNext[ip], xp[3*ip], xp[3*ip+1], xp[3*ip+2], distP, rCoordPPrevious[ip]);
        }
    }else if(spi_velocity_mode==OptionConstants::EQTERM_SPI_VELOCITY_MODE_NONE){
        for(len_t ip=0;ip<nShard;ip++){
        	// If the shards do not move, we can not use the distance 
        	// the shards travel in one time step as a length scale to set the tolerance.
        	// Here we use a hardcoded length scale of 1 cm
            if(t!=t_old)
                rGrid->GetRThetaPhiFromCartesian(&rCoordPPrevious[ip], &thetaCoordPPrevious[ip], &phiCoordPPrevious[ip], xpPrevious[3*ip], xpPrevious[3*ip+1], xpPrevious[3*ip+2], 0.01, rCoordPPrevious[ip]);
                
            rGrid->GetRThetaPhiFromCartesian(&rCoordPNext[ip], &thetaCoordPNext[ip], &phiCoordPNext[ip], xp[3*ip], xp[3*ip+1], xp[3*ip+2], 0.01, rCoordPPrevious[ip]);
        }
    }else {throw DREAMException("SPIHandler: unrecognized SPI shard velocity mode");}
    
    // Calculate the radial index of each shard
    CalculateIrp();
    
    // Calculate ablation rate (if any)
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        CalculateYpdotNGSParksTSDW();
        if(t!=t_old)
            for(len_t ip=0; ip<nShard; ip++)
                YpdotPrevious[ip] = Ypdot[ip];       
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){ 
    
        for(len_t ir=0;ir<nr;ir++){
            // Total electron heat flux. 
            // The factor 1/4 is an approximate way to convert from the flux in all directions to the flux in only one direction
            qtot[ir]=(qhot[ir] + 4.0*sqrt(2.0/(M_PI*Constants::me))*ncold[ir]*pow(Constants::ec*Tcold[ir],3.0/2.0))/4.0;
            
            // Effective energy of incomming electrons
            Eeff[ir]=4.0/3.0*(Wcold[ir]+Whot[ir])/ntot[ir];
        }
        CalculateYpdotNGSParksTSDWKinetic();
        if(t!=t_old)
            for(len_t ip=0; ip<nShard; ip++)
                YpdotPrevious[ip] = Ypdot[ip]; 
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT){
        for(len_t ip=0;ip<nShard;ip++){
            Ypdot[ip]=0;
            YpdotPrevious[ip]=0;
        }
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_NGPS){
        throw NotImplementedException("SPIHandler: NGPS ablation is not yet supported");
    }else {throw DREAMException("SPIHandler: unrecognized SPI shard ablation mode");}
    
    // Calculate magnetic field damping (if any)
    if(spi_magnetic_field_dependence_mode==OptionConstants::EQTERM_SPI_MAGNETIC_FIELD_DEPENDENCE_MODE_JOREK)
        for(len_t ip = 0; ip<nShard; ip++){
            Ypdot[ip]*=CalculateBFieldDampingJOREK(irp[ip]);
            if(t!=t_old)
                YpdotPrevious[ip]*=CalculateBFieldDampingJOREK(irp[ip]);
        }
    // Calculate radius of the neutral cloud (if any)
    if(spi_cloud_radius_mode!=OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_NEGLECT)
        CalculateRCld();

    // Calculate deposition (if any)
    if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL || spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){
        CalculateTimeAveragedDeltaSourceLocal(depositionProfilesAllShards);
        
    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        CalculateGaussianSourceLocal(depositionProfilesAllShards);

    }else if(spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT){
        for(len_t ir=0;ir<nr;ir++)
            depositionRate[ir]=0;
    }else {throw DREAMException("SPIHandler: unrecognized SPI material deposition mode");}
    
    // Calculate drift (if any)
    if(spi_shift_mode!=OptionConstants::EQTERM_SPI_SHIFT_MODE_NEGLECT || spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){
    
        // Calculate the projection of the major radius unit vector on the flux surface normal, cos(thetaDrift)=R dot grad r / (|R|*|grad r|).
        // The drift is in general not only directed along the major radius, as the E-field causing the drift rotates
        // with the twist of the field lines, but it turns out that the radial drift component will anyway be modified 
        // by this factor. This can be seen by noting that the angle between Yhat and yhat in Vallhagen et al JPP 2023
        // will everywhere be modified by thetaDrift for shards not on the outboard midplane, and after performing the 
        // integral in eq. 2.8 an overall modifying factor of cos(thetaDrift) falls out.
        if(t!=t_old){
            for(len_t ip=0;ip<nShard;ip++){
                if(rCoordPPrevious[ip]<rGrid->GetR_f(nr)){
                    rGrid->GetGradRCartesian(gradRCartesianPrevious,rCoordPPrevious[ip],thetaCoordPPrevious[ip],phiCoordPPrevious[ip]);
                    
                    // Note that even if we need a finite value for the major radius when calculating the drift, this value
                    // does not affect the angle against the flux surface normal, if the major radius of the radial grid is infinite.
                    // In that case, the major radius unit vector simply coincides with the unit vector for x in the SPI coordinate
                    // system.
                    real_t R0 = this->rGrid->GetR0();
                    if(isinf(R0))
                        cosThetaDrift[ip] = gradRCartesianPrevious[0]/sqrt(gradRCartesianPrevious[0]*gradRCartesianPrevious[0] + gradRCartesianPrevious[1]*gradRCartesianPrevious[1] + gradRCartesianPrevious[2]*gradRCartesianPrevious[2]);
                    else
                        cosThetaDrift[ip] = ((xpPrevious[3*ip] + R0)*gradRCartesianPrevious[0] + xpPrevious[3*ip+2]*gradRCartesianPrevious[2])/sqrt(((xpPrevious[3*ip] + R0)*(xpPrevious[3*ip] + R0) + xpPrevious[3*ip+2]*xpPrevious[3*ip+2]) * (gradRCartesianPrevious[0]*gradRCartesianPrevious[0] + gradRCartesianPrevious[1]*gradRCartesianPrevious[1] + gradRCartesianPrevious[2]*gradRCartesianPrevious[2]));
                } else {
                    cosThetaDrift[ip]=1;
                }
            }
        }
    
        if(spi_shift_mode==OptionConstants::EQTERM_SPI_SHIFT_MODE_ANALYTICAL){
            // We only calculate the drift once per time step, to avoid a discontinuity between Newton iterations
            if(t!=t_old){
                for(len_t ip=0;ip<nShard;ip++){
                    if (YpPrevious[ip]>0 && irp[ip]<nr){
                        YpConversion(ip);
                        AssignShardSpecificDriftParameters(ip);
                        AssignDriftComputationParameters(ip);
                        AssignDriftTimeParameters(ip);
                        shift_r[ip] = Deltar(ip);
                        nbrShiftGridCell[ip] = CalculateDriftIrp(ip, shift_r[ip]);// Negative if shift is towards smaller radii
                        shift_store[ip] = shift_r[ip];
                    }else{
                        nbrShiftGridCell[ip]=0;
                        shift_store[ip]=0;
                    }
                }
            }
        } else if(spi_shift_mode==OptionConstants::EQTERM_SPI_SHIFT_MODE_PRESCRIBED || spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){// Prescribed drift (in terms of grid cells)
            // If the radial projection of the drift is negative, the shift should go towards smaller radii, 
            // unless the shift goes past the core and ends at a larger radii on the other side. Here we 
            // account for this when setting the prescribed shift, keeping in mind that nbrShiftGridCell 
            // should be negative if the shift is towards smaller radii, while nbrShiftGridCellPrescribed 
            // is always positive
            for(len_t ip=0;ip<nShard;ip++){
                if(cosThetaDrift[ip]<0){// Is the radial projection of the drift negative?
                    nbrShiftGridCell[ip] = std::abs((int_t)irp[ip]-nbrShiftGridCellPrescribed[ip])-(int_t)irp[ip];
                    
                    //Account for that grid cell 0 should be counted twice (on both sides of the magnetic axis)
                    if((int_t)irp[ip]-nbrShiftGridCellPrescribed[ip]<0)
                        nbrShiftGridCell[ip]--;    
                }else if(rCoordPNext[ip]<rCoordPPrevious[ip]){
                    nbrShiftGridCell[ip] = nbrShiftGridCellPrescribed[ip];
                }
            }
        }
        // Shift the deposition profile
        for(len_t ip=0;ip<nShard;ip++){
            if(nbrShiftGridCell[ip]<0){
                for(len_t ir=0;ir<nr;ir++){
                    if(ir>=nr+nbrShiftGridCell[ip])
                        depositionProfilesAllShards[ir*nShard+ip]=0;
                    else
                        // recall that nbrShiftGridCell is now negative here, 
                        // so the point that we move inwards, with index ir-nbrShiftGridCell, is really outside the point with index ir
                        depositionProfilesAllShards[ir*nShard+ip]=rGrid->GetVpVol((int_t)ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*depositionProfilesAllShards[((int_t)ir-nbrShiftGridCell[ip])*nShard+ip];
                }
            }else if(nbrShiftGridCell[ip]>0){
                for(int_t ir=nr-1;ir>=0;ir--){// Use ant int as a loop index so that the loop can be terminated by ir becoming <0
                    if(ir<nbrShiftGridCell[ip])
                        depositionProfilesAllShards[ir*nShard+ip]=0;
                    else
                        depositionProfilesAllShards[ir*nShard+ip]=rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*depositionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                }
            }
        }
    }else if(spi_shift_mode==OptionConstants::EQTERM_SPI_SHIFT_MODE_NEGLECT){
        for(len_t ip=0;ip<nShard;ip++){
            nbrShiftGridCell[ip] = 0;
        }
    }else {throw DREAMException("SPIHandler: unrecognized SPI shift mode");}

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
    this->t_old=t;	
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


real_t SPIHandler::CalculateBFieldDampingJOREK(len_t ir){
    if(ir<rGrid->GetNr()){
        if((rGrid->GetFSA_B(ir)*rGrid->GetBmin(ir))>2.0)
            return pow(2.0/(rGrid->GetFSA_B(ir)*rGrid->GetBmin(ir)),0.843);
        else
            return 1.0;
    }else{
        return 0.0;
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
            if(YpPrevious[ip]>0 && irp[ip]<nr){
                real_t heatAbsorbtionPrefactor = plasmoidAbsorbtionFactor[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]/exp(1.8);// Dividing by exp(1.8) as an approximate way to account for the sheath potential damping of the absorbed heat flux (see e.g. Parks et al PoP 2000, https://doi.org/10.1063/1.874052)
                
                heatAbsorbtionRate[ir]+=-heatAbsorbtionPrefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                
                // Account for shifted re-deposition 
                // NOTE: only strictly valid for delta function kernel (assumes deposition only on one side of r=0)
                if(nbrShiftGridCell[ip]<0){
                    if(ir<nr+nbrShiftGridCell[ip])
                        heatAbsorbtionRate[ir]+=rGrid->GetVpVol((int_t)ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*heatAbsorbtionPrefactor*heatAbsorbtionProfilesAllShards[((int_t)ir-nbrShiftGridCell[ip])*nShard+ip];
                }else if(nbrShiftGridCell[ip]>=0){
                    if((int_t)ir>=nbrShiftGridCell[ip])
                        heatAbsorbtionRate[ir]+=rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*heatAbsorbtionPrefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                }     
            }
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

        if(irp[ip] < nr){
            // Find out if the shard has passed its turning point (where the radial coordinate goes from decreasing to increasing)
            // If so, we split the averaging in two steps, before and after the turning point (where the analythical expression used for the averaging breaks down)
            // We determine wether the turning point has been passed by checking if there is a sign change of the 
            // projection of the direction vector for the shard motion on the gradient of the radial coordinate before and after the time step
            nSplit=1;
            iSplit=0;
            turningPointPassed=false;
            rGrid->GetGradRCartesian(gradRCartesian,rCoordPNext[ip],thetaCoordPNext[ip],phiCoordPNext[ip]);
            rGrid->GetGradRCartesian(gradRCartesianPrevious,rCoordPPrevious[ip],thetaCoordPPrevious[ip],phiCoordPPrevious[ip]);
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
		        if(rSourceMax!=rSourceMin){
                	for(len_t ir=0;ir<nr;ir++){
      
                    	if(!(rGrid->GetR_f(ir)>rSourceMax || rGrid->GetR_f(ir+1)<rSourceMin)){
		                    timeAveragedDeltaSource[ir*nShard+ip]+=1.0/(rGrid->GetVpVol(ir)*VpVolNormFactor*(rSourceMax-rSourceMin))*
		                                                         (min(rGrid->GetR_f(ir+1),rSourceMax)-max(rGrid->GetR_f(ir),rSourceMin))/rGrid->GetDr(ir); 
                                              
		                }
		            }
                }
                iSplit++;
            }
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
* Function used to calculate the number of grid cells to shift the deposition due to the drift
* (the shift is only made in integer steps of the radial resolution)
* ip: shard index
* shift: the radial shift to be made (exact, not necessarilly an integer of radial grid cells), sign included
*/
int_t SPIHandler::CalculateDriftIrp(len_t ip, real_t shift){
    for(len_t ir=0; ir<nr;ir++){
        if(abs(rCoordPNext[ip] + shift)<rGrid->GetR_f(ir+1) && abs(rCoordPNext[ip] + shift)>rGrid->GetR_f(ir)){
            return (int_t)ir - (int_t)irp[ip];
        }
    }
    return nr;
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
bool SPIHandler::setJacobianYpdot(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_FLUID_NGS){
        return setJacobianYpdotNGS(jac, derivId, scaleFactor);
    }else if(spi_ablation_mode==OptionConstants::EQTERM_SPI_ABLATION_MODE_KINETIC_NGS){
        return setJacobianYpdotNGSKinetic(jac, derivId, scaleFactor);
    }else{
        return false;
    }
}

/**
 * Wrapper function for calculation of partial derivatives of the deposition rate, depending on the settings
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS and to apply weights (eg ionization energies and/or molar fractions)
 */
bool SPIHandler::setJacobianDepositionRate(FVM::Matrix *jac, len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset){
    if((spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL || 
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE)||
        spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN){
        return setJacobianDepositionRateDensCons(jac, derivId, scaleFactor, SPIMolarFraction, rOffset);
    }else{
        return false;
    }
}

/**
 * Wrapper function for calculation of partial derivatives of the heat absorbtion rate, depending on the settings
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
bool SPIHandler::setJacobianAdiabaticHeatAbsorbtionRate(FVM::Matrix *jac, len_t derivId, real_t scaleFactor){
    if(spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS ||
        spi_heat_absorbtion_mode==OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN){
        return setJacobianAdiabaticHeatAbsorbtionRateMaxwellian(jac, derivId, scaleFactor);
    }else{
        return false;
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
bool SPIHandler::setJacobianYpdotNGS(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    bool jacIsSet=false;
    if(derivId==id_Tcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*5.0/3.0*Ypdot[ip]/Tcold[irp[ip]]);
                jacIsSet=true;
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/ncold[irp[ip]]);
                jacIsSet=true;
            }
        }
    }
    return jacIsSet;
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
bool SPIHandler::setJacobianYpdotNGSKinetic(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    bool jacIsSet=false;
    if(derivId==id_Tcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){        		
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]*3.0/2.0*(qtot[irp[ip]]-qhot[irp[ip]])/Tcold[irp[ip]]);
                jacIsSet=true;
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]*(qtot[irp[ip]]-qhot[irp[ip]])/ncold[irp[ip]]);
                jacIsSet=true;
            }
        }
    }else if(derivId==id_Whot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*7.0/6.0*Ypdot[ip]/(Whot[irp[ip]]+Wcold[irp[ip]]));
                jacIsSet=true;
            }
        }
    }else if(derivId==id_qhot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*1.0/3.0*Ypdot[ip]/qtot[irp[ip]]);
                jacIsSet=true;
            }
        }
    }else if(derivId==id_ntot){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], -scaleFactor*7.0/6.0*Ypdot[ip]/ntot[irp[ip]]);
                jacIsSet=true;
            }
        }
    }else if(derivId==id_Wcold){
        for(len_t ip=0;ip<nShard;ip++){
            if(irp[ip]<nr){
                jac->SetElement(ip,irp[ip], scaleFactor*7.0/6.0*Ypdot[ip]/(Whot[irp[ip]]+Wcold[irp[ip]]));
                jacIsSet=true;
            }
        }
    }
    return jacIsSet;
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
bool SPIHandler::setJacobianDepositionRateDensCons(FVM::Matrix *jac,len_t derivId, real_t *scaleFactor, real_t *SPIMolarFraction, len_t rOffset){
    bool jacIsSet=false;
    if(derivId==id_Yp){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(YpPrevious[ip]>0){
                    jac->SetElement(ir+rOffset,ip,-scaleFactor[ir]*SPIMolarFraction[ip]*12.0/5.0*M_PI*pow(abs(Yp[ip]),4.0/5.0)/pelletMolarVolume[ip]*Constants::N_Avogadro/dt*depositionProfilesAllShards[ir*nShard+ip]);
                    jacIsSet=true;
                }
            }
        }
    }
    return jacIsSet;
}

/**
 * Sets jacobian block contribution from the heat absorbtion rate,
 * assuming Maxwellian distribution for the insident electrons
 *
 * jac: jacobian block to set partial derivatives to
 * derivId: ID for variable to differentiate with respect to
 * scaleFactor: Used to move terms between LHS and RHS
 */
bool SPIHandler::setJacobianAdiabaticHeatAbsorbtionRateMaxwellian(FVM::Matrix *jac,len_t derivId, real_t scaleFactor){
    bool jacIsSet=false;
    if(derivId==id_Yp){
        if(spi_cloud_radius_mode==OptionConstants::EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT){
            for(len_t ir=0;ir<nr;ir++){
                for(len_t ip=0;ip<nShard;ip++){
                    if(YpPrevious[ip]>0 && irp[ip]<nr){
                        real_t prefactor = -scaleFactor*6.0/5.0/Yp[ip]*plasmoidAbsorbtionFactor[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]/exp(1.8);
                        real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                        // Account for shifted re-deposition
                        if(nbrShiftGridCell[ip]<0 && ir<nr+nbrShiftGridCell[ip])
                            jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                        else if(nbrShiftGridCell[ip]>0 && ir>=(len_t)nbrShiftGridCell[ip])
                            jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                    
                        jac->SetElement(ir,ip,jacEl);
                        jacIsSet=true;
                    }
                }
            }
        }
    }else if(derivId==id_Tcold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(irp[ip]<nr){
                    real_t prefactor = -scaleFactor*3.0/2.0*plasmoidAbsorbtionFactor[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec/exp(1.8);
                    real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                    // Account for shifted re-deposition
                    if(nbrShiftGridCell[ip]<0 && ir<nr+nbrShiftGridCell[ip])
                        jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                    else if(nbrShiftGridCell[ip]>0 && ir>=(len_t)nbrShiftGridCell[ip])
                        jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                
                    jac->SetElement(ir,irp[ip],jacEl);
                    jacIsSet=true;
                }
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(irp[ip]<nr){
                    real_t prefactor = -scaleFactor*plasmoidAbsorbtionFactor[ip]*M_PI*rCld[ip]*rCld[ip]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]/exp(1.8);
                    real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                    // Account for shifted re-deposition
                    if(nbrShiftGridCell[ip]<0 && ir<nr+nbrShiftGridCell[ip])
                        jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                    else if(nbrShiftGridCell[ip]>0 && ir>=(len_t)nbrShiftGridCell[ip])
                        jacEl+=-rGrid->GetVpVol(ir-nbrShiftGridCell[ip])/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-nbrShiftGridCell[ip])*nShard+ip];
                
                    jac->SetElement(ir,irp[ip],jacEl);
                    jacIsSet=true;
                }
            }
        }
    }
    return jacIsSet;
}
