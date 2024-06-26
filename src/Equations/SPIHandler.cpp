/**
 * Implementation of a class that calculates and stores quantities related to the SPI shards
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include "DREAM/Equations/SPIHandler.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>

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
const real_t Zavg0 = 1;// Assume all ion species are singly ionized directly after the neutral phase and while the cloud detaches from the pellet
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
    real_t *T, real_t T_0, real_t delta_y, real_t Rm, real_t *ZavgArray, len_t nZavg, real_t *Zs, real_t *isotopesDrift,
    real_t VpVolNormFactor=1, real_t rclPrescribedConstant=0.01, int_t *nbrShiftGridCell=nullptr){

    // Get pointers to relevant objects
    this->rGrid=g->GetRadialGrid();
    this->unknowns=u;
    this->VpVolNormFactor=VpVolNormFactor;

	// Get the major radius, to be used to properly normalize VpVol
	real_t R0 = this->rGrid->GetR0();
    rf = this->rf;
    q = 1; //TODO

	// If R0 is infinite, i.e. toroidicity is not included in the simulation,
	// we can not use R0 from the radial grid of this simulation to calculate 
	// the size of the flux surfaces. The corresponding factor correcting the 
	// size of the flux surfaces must instead be included directly in the 
	// VpVolNormFactor
	if(!isinf(R0)){
	    this->VpVolNormFactor*=R0;
        if(Rm==-1)
            Rm=R0;
    }else if(Rm==-1 && spi_shift_mode==true)
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

    // Memory allocation
    AllocateQuantities();
    
    this->T=T;
    this->T_0=T_0;
    this->delta_y=delta_y;
    this->Rm=Rm;
    this->NZ=NZ;
    this->ZavgArray=ZavgArray;
    
    // Initialize rCoordPrevious to the radial coordinate at the plasma edge 
    // to use as a starting guess if rCoord must be solved for numerically
    for(len_t ip=0;ip<nShard;ip++)
        rCoordPPrevious[ip]=rGrid->GetR_f(nr-1);

    // Calculate pellet molar mass, molar volume and density
    real_t molarMass=0;
    real_t solidDensity=0;
    real_t ZavgList=0;
    for(len_t ip=0;ip<nShard;ip++){
        pelletMolarMass[ip]=0;
        pelletMolarVolume[ip]=0;
        Zavg[ip]=0;
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
            for(len_t i=0;i<nZavg;i++){
                if(Z[iZ]==Zs[i] && isotopes[iZ]==isotopesDrift[i]){
                    ZavgList=ZavgArray[i];
                }
            }
            for(len_t ip=0;ip<nShard;ip++){
		        pelletMolarMass[ip]+=molarMass*molarFraction[offset+ip];
		        pelletMolarVolume[ip]+=molarMass/solidDensity*molarFraction[offset+ip];
                
		        if(Z[iZ]==1 && isotopes[iZ]==2)
		            pelletDeuteriumFraction[ip]+=molarFraction[offset+ip];

                Zavg[ip]+=ZavgList*molarFraction[offset+ip];
                if (!((Z[iZ]==1 && isotopes[iZ]==2) || Z[iZ]==10) && counter2==0){
                    printf("SPIHandler: Using other molarFractions other than deuterium and neon will lead to urealistic ablating.");
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
	else if(nbrShiftGridCell!=nullptr)
	    for(len_t ip=0;ip<nShard;ip++){
	        this->nbrShiftGridCellPrescribed[ip] = nbrShiftGridCell[ip];
        }
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
    T = new real_t[nShard];
    pelletDeuteriumFraction=new real_t[nShard];
    rp=new real_t[nShard];
    rpdot=new real_t[nShard];
    shift_r=new real_t[nShard];
    shift_store=new real_t[nShard];
    YpdotPrevious=new real_t[nShard];
    Zavg=new real_t[nShard];
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
    delete [] T;
    delete [] pelletDeuteriumFraction;
    delete [] rp;
    delete [] rpdot;
    delete [] shift_r;
    delete [] shift_store;
    delete [] YpdotPrevious;
    delete [] Zavg;
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
void SPIHandler::AssignShardSpecificParameters(len_t ip){
    v0 = -vp[3*ip];
    n_e = ncoldPrevious[irp[ip]];
    Te = TcoldPrevious[irp[ip]];
    B = sqrt(this->rGrid->GetFSA_B2(irp[ip])) * rGrid->GetBmin(irp[ip]);
    sigma = rf->GetElectricConductivity(irp[ip]);
    n_i = 0;
    for(len_t iZ=0;iZ<NZ;iZ++){
        n_i += rf->GetIonHandler()->GetTotalIonDensity(irp[ip], iZ);
    }
}

/**
 * Computes relevant quantities for the drift displacement of a shard
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
void SPIHandler::AssignComputationParameters(len_t ip){
    CST = sqrt((gamma_e*Zavg[ip] + gamma_i) * qe * T[ip]/(pelletMolarMass[ip]/N_Avogadro));
    CST0 = sqrt((gamma_e*Zavg[ip] + gamma_i) * qe * T_0/(pelletMolarMass[ip]/N_Avogadro));
    G = -4 * M_PI * pelletDensity[ip] * rp[ip] * rp[ip] * rpdot[ip];
    n_0 = (1 + Zavg0)*G/(2 * M_PI * delta_y * delta_y * pelletMolarMass[ip]/N_Avogadro * CST0);
    a0 = ((1 + Zavg0)*qe*T_0/(pelletMolarMass[ip]/N_Avogadro*Rm));
    t_detach = -v0/a0 + sqrt(v0*v0/(a0*a0) + 2*delta_y/a0);
    Lc = 2 * CST0*t_detach;
    n = n_0 * Lc; 
    v_lab = a0 * t_detach;
    Reff = -2*M_PI*M_PI*Rm*this->rGrid->GetMinorRadius()/(sigma*delta_y*delta_y*delta_y*log(delta_y/(this->rGrid->GetMinorRadius()*M_PI)));
}

/**
 * Computes the characteristic quantities of time for a shard
 * t_acc : Acceleration time scale
 * t_pol : Poloidal rotation time scale
 * t_pe  : Pressure equilibration time scale
 * t_exp : Expansion time scale
 * t_polp: Normalized t_pol
 * t_pep : Normalized t_pe
 * t_expp: Normalized t_exp
 */
void SPIHandler::AssignTimeParameters(len_t ip){
    t_acc = n/(1+Zavg[ip])*pelletMolarMass[ip]/N_Avogadro*Reff/(B*B);
    t_pol = q*Rm/CST;
    t_pe = qe*T[ip]*n/(2*CST*(n_i+n_e)*qe*Te);
    t_exp = Lc/(2*CST);
    t_polp = t_pol/t_acc;
    t_pep = t_pe/t_acc;
    t_expp = t_exp/t_acc;
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


real_t SPIHandler::BisFunction(real_t t_prim){
    return t_prim + t_expp;
}

//Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveFirstRow(real_t t_prim){//- at cos in article
    real_t t_bis = BisFunction(t_prim);
    real_t term1 = Epsiloni(t_bis, t_bis/t_polp);
    real_t numerator = (sin(t_bis/t_polp) + 1/t_polp * cos(t_bis/t_polp));
    real_t denominator =  t_pep * (1 + 1/(t_polp * t_polp));
    real_t term2 = numerator/denominator;
    return term1 - term2;
}

// Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveSecondRow(real_t t_prim){
    real_t t_bis = BisFunction(t_prim);
    real_t term1 = Epsiloni(0, t_bis/t_polp);
    real_t term2 = Epsiloni(t_bis, t_bis/t_polp);
    real_t result = term1 - term2;
    return result;
}

// Function to evaluate two of the terms in equation (A4) 
real_t SPIHandler::PrimitiveThirdRow(real_t t_prim){
    real_t t_bis = BisFunction(t_prim);
    real_t term1 = t_polp * cos(t_bis/t_polp);
    real_t term2 = sin(t_bis/t_polp);
    return term1 + term2;
}

real_t SPIHandler::FirstRow(){
    real_t term1 = PrimitiveFirstRow(t_pep);
    real_t term2 = PrimitiveFirstRow(0);
    return term1 - term2;
}

real_t SPIHandler::SecondRow(){
    real_t term1 = PrimitiveSecondRow(t_pep);
    real_t term2 = PrimitiveSecondRow(0);
    return term1 - term2;
}

real_t SPIHandler::ThirdRow(){
    real_t factor = 1/t_pep * 1/(1 + 1/(t_polp*t_polp));
    real_t term1 = PrimitiveThirdRow(t_pep);
    real_t term2 = PrimitiveThirdRow(0);
    return factor * (term1 - term2);
}

// Function to collect all terms to evaluate equation A4
real_t SPIHandler::Deltar(len_t ip){
    real_t first = FirstRow();
    real_t second = SecondRow();
    real_t third = ThirdRow();
    real_t term1 = v_lab * t_acc;
    real_t factor = (1+Zavg[ip])*2*qe*T[ip]*q/(CST*pelletMolarMass[ip]/N_Avogadro)*t_acc;
    return (term1 + factor * (first + second + third))*cos(thetaCoordPPrevious[ip]);
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
        // Calculate drift (if any)
        if(spi_shift_mode==OptionConstants::EQTERM_SPI_SHIFT_MODE_ANALYTICAL){
            if(t!=t_old){
                for(len_t ip=0;ip<nShard;ip++){
                    if (YpPrevious[ip]>0 && irp[ip]<nr){
                        YpConversion(ip);
                        AssignShardSpecificParameters(ip);
                        AssignComputationParameters(ip);
                        AssignTimeParameters(ip);
                        shift_r[ip] = Deltar(ip);
                        nbrShiftGridCell[ip] = CalculateDriftIrp(ip, shift_r[ip]);
                        shift_store[ip] = shift_r[ip];
                    }
                }
            }
        } else if(nbrShiftGridCellPrescribed!=nullptr || spi_deposition_mode==OptionConstants::EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE){// Prescribed drift (in terms of grid cells)
            for(len_t ip=0;ip<nShard;ip++){
                if(rCoordPNext[ip]>rCoordPPrevious[ip]){
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
                real_t heatAbsorbtionPrefactor = M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]]/exp(1.8);// Dividing by exp(1.8) as an approximate way to account for the sheath potential damping of the absorbed heat flux
                
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

int_t SPIHandler::CalculateDriftIrp(len_t ip, real_t shift){
    for(len_t ir=0; ir<nr;ir++){
        if(rCoordPNext[ip]>rCoordPPrevious[ip] && abs(rCoordPNext[ip] - shift)<rGrid->GetR_f(ir+1) && abs(rCoordPNext[ip] - shift)>rGrid->GetR_f(ir)){
            return (int_t)ir - (int_t)irp[ip];
        }else if(rCoordPNext[ip] + shift<rGrid->GetR_f(ir+1) && rCoordPNext[ip] + shift>rGrid->GetR_f(ir)){
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
                        real_t prefactor = -scaleFactor*6.0/5.0/Yp[ip]*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]];
                        real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                        // Account for shifted re-deposition
                        if(rCoordPNext[ip]>rCoordPPrevious[ip] && ir<nr-1)
                            jacEl+=-rGrid->GetVpVol(ir+1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir+1)*nShard+ip];
                        else if(rCoordPNext[ip]<rCoordPPrevious[ip] && ir>0)
                            jacEl+=-rGrid->GetVpVol(ir-1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-1)*nShard+ip];
                    
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
                    real_t prefactor = -scaleFactor*3.0/2.0*M_PI*rCld[ip]*rCld[ip]*ncold[irp[ip]]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec;
                    real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                    // Account for shifted re-deposition
                    if(rCoordPNext[ip]>rCoordPPrevious[ip] && ir<nr-1)
                        jacEl+=-rGrid->GetVpVol(ir+1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir+1)*nShard+ip];
                    else if(rCoordPNext[ip]<rCoordPPrevious[ip] && ir>0)
                        jacEl+=-rGrid->GetVpVol(ir-1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-1)*nShard+ip];
                
                    jac->SetElement(ir,irp[ip],jacEl);
                    jacIsSet=true;
                }
            }
        }
    }else if(derivId==id_ncold){
        for(len_t ir=0;ir<nr;ir++){
            for(len_t ip=0;ip<nShard;ip++){
                if(irp[ip]<nr){
                    real_t prefactor = -scaleFactor*M_PI*rCld[ip]*rCld[ip]*sqrt(8.0*Constants::ec*Tcold[irp[ip]]/(M_PI*Constants::me))*Constants::ec*Tcold[irp[ip]];
                    real_t jacEl = prefactor*heatAbsorbtionProfilesAllShards[ir*nShard+ip];
                        
                    // Account for shifted re-deposition
                    if(rCoordPNext[ip]>rCoordPPrevious[ip] && ir<nr-1)
                        jacEl+=-rGrid->GetVpVol(ir+1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir+1)*nShard+ip];
                    else if(rCoordPNext[ip]<rCoordPPrevious[ip] && ir>0)
                        jacEl+=-rGrid->GetVpVol(ir-1)/rGrid->GetVpVol(ir)*prefactor*heatAbsorbtionProfilesAllShards[(ir-1)*nShard+ip];
                
                    jac->SetElement(ir,irp[ip],jacEl);
                    jacIsSet=true;
                }
            }
        }
    }
    return jacIsSet;
}
