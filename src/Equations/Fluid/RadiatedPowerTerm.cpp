#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"


/**
 * Implementation of a class which represents the 
 * radiated power as calculated with rate coefficients
 * from the ADAS database (PLT corresponds to line
 * and PRB to brems and recombination radiated power).
 * The term is of the form n_e * sum_i n_i L_i, summed over all
 * ion species i. In the semi-implicit solver, n_e is the "unknown"
 * evaluated at the next time step and n_i L_i coefficients.
 * We ignore the Jacobian with respect to L_i(n,T) and capture only the
 * n_e and n_i contributions.
 */


using namespace DREAM;


RadiatedPowerTerm::RadiatedPowerTerm(
    FVM::Grid* g, FVM::UnknownQuantityHandler *u, IonHandler *ionHandler, 
	ADAS *adas, NIST *nist, AMJUEL* amjuel,
    enum OptionConstants::ion_opacity_mode *opacity_modes, bool includePRB
) : FVM::DiagonalComplexTerm(g,u), includePRB(includePRB) {

    SetName("RadiatedPowerTerm");

    this->adas = adas;
    this->nist = nist;
    this->amjuel = amjuel;
    this->ionHandler = ionHandler;
    
    this->opacity_modes = new enum OptionConstants::ion_opacity_mode[ionHandler->GetNZ()];
    for(len_t iz=0;iz<ionHandler->GetNZ();iz++)
    	this->opacity_modes[iz] = opacity_modes[iz];

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    AddUnknownForJacobian(unknowns, id_ncold);
    AddUnknownForJacobian(unknowns, id_ni);
    AddUnknownForJacobian(unknowns, id_Tcold);

    // constants appearing in bremsstrahlung formula
    real_t c = Constants::c;
    this->bremsPrefactor = (32.0/3.0)*Constants::alpha*Constants::r0*Constants::r0*c
        * sqrt(Constants::me*c*c*Constants::ec*2.0/M_PI);
    this->bremsRel1 = 19.0/24.0; // relativistic-maxwellian correction
    this->bremsRel2 = 5.0/(8.0*M_SQRT2)*(44.0-3.0*M_PI*M_PI); // e-e brems correction
}


/**
 * Set the weights of this term.
 */
void RadiatedPowerTerm::SetWeights(){
    this->SetWeights(1.0);
}

/**
 * Set the weights of this term, scaling contributions from neutral particles
 * by the given factor 'neutralScale'.
 *
 * neutralScale: Scale factor to multiply contributions from neutral particles with.
 */
void RadiatedPowerTerm::SetWeights(const real_t neutralScale) {
    len_t NCells = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();
    
    real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t *n_i    = unknowns->GetUnknownData(id_ni);
    
    for (len_t i = 0; i < NCells; i++)
            weights[i] = 0;

    for(len_t iz = 0; iz<nZ; iz++){
        ADASRateInterpolator *PLT_interper = adas->GetPLT(Zs[iz]);
        ADASRateInterpolator *PRB_interper = adas->GetPRB(Zs[iz]);
        ADASRateInterpolator *ACD_interper = adas->GetACD(Zs[iz]);
        ADASRateInterpolator *SCD_interper = adas->GetSCD(Zs[iz]);
        real_t dWi = 0;
        real_t Li = 0;
        real_t Bi = 0;
        for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
            len_t indZ = ionHandler->GetIndex(iz,Z0);
            for (len_t i = 0; i < NCells; i++){

            	if(Zs[iz]==1 && opacity_modes[iz]==OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){//Ly-opaque deuterium radiation from AMJUEL
		            // Radiated power term
		            Li = amjuel->getIonizLossLyOpaque(Z0, n_cold[i], T_cold[i]);// includes both line radiation and ionization potential energy difference
		            
		            // The AMJUEL coefficients for recombination radiation do not contain bremsstrahlung,
		            // but are on the other hand adjusted for repeated excitation/deexcitation and three-body recombination
		            // Thus, the recombination radiation and recombination gain binding energy term should be included
		            // regardless of wether includePRB is true or false
		            Li += amjuel->getRecRadLyOpaque(Z0, n_cold[i], T_cold[i]);
		            
		            Bi = 0;
		            // Binding energy rate term
		            if(Z0>0){       // Recombination gain
		            	dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0-1);
		                Bi -= dWi * amjuel->getRecLyOpaque(Z0, n_cold[i], T_cold[i]);
		                
		                // As the AMJUEL coefficients do not include bremsstrahlung,
		                // if this is expected to be a part of Li we need to add it explicitly.
		                // We do this in a similar way as would have been done for all species if includePRB==false
    	                if(includePRB)
	                        Li+=bremsPrefactor*sqrt(T_cold[i])*Z0*Z0*(1 + bremsRel1*T_cold[i]/Constants::mc2inEV);
	                }
            	}else{
		            // Radiated power term
		            Li =  PLT_interper->Eval(Z0, n_cold[i], T_cold[i]);
		            if (includePRB) 
		                Li += PRB_interper->Eval(Z0, n_cold[i], T_cold[i]);
		            Bi = 0;
		            // Binding energy rate term
		            if(Z0>0 && includePRB)       // Recombination gain
		                Bi -= dWi * ACD_interper->Eval(Z0, n_cold[i], T_cold[i]);
		            if(Z0<Zs[iz]){ // Ionization loss
		                dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0);
		                Bi += dWi * SCD_interper->Eval(Z0, n_cold[i], T_cold[i]);
		            }

                }

                if (Z0 == 0)
                    weights[i] += neutralScale * n_i[indZ*NCells + i]*(Li+Bi);
                else
                    weights[i] += n_i[indZ*NCells + i]*(Li+Bi);
            }
        }
    }

    /**
     * If neglecting the recombination radiation, explicitly add the
     * bremsstrahlung loss which is otherwise included in 'PRB'.
     * Using the R J Gould, The Astrophysical Journal 238 (1980)
     * formula with a relativistic correction. The 'bremsPrefactor'
     * is the factor ~1.69e-38 appearing in the NRL formula (implemented in GO).
     * The 'ionTerm' equals sum_ij n_i^(j) Z_0j^2 
     */
    if(!includePRB) 
        for(len_t i=0; i<NCells; i++){
            real_t ionTerm = ionHandler->GetZeff(i)*ionHandler->GetFreeElectronDensityFromQuasiNeutrality(i);
            real_t relativisticCorrection = bremsRel1*ionTerm + bremsRel2*n_cold[i];
            weights[i] += bremsPrefactor*sqrt(T_cold[i])*(ionTerm + relativisticCorrection*T_cold[i]/Constants::mc2inEV);
        }
}

void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t indZs){
    this->SetDiffWeights(derivId, indZs, 1.0);
}

void RadiatedPowerTerm::SetDiffWeights(len_t derivId, len_t, const real_t neutralScale) {
    len_t NCells = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();

    real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t *n_i    = unknowns->GetUnknownData(id_ni);

    if(derivId == id_ni){
        for(len_t iz = 0; iz<nZ; iz++){
            ADASRateInterpolator *PLT_interper = adas->GetPLT(Zs[iz]);
            ADASRateInterpolator *PRB_interper = adas->GetPRB(Zs[iz]);
            ADASRateInterpolator *ACD_interper = adas->GetACD(Zs[iz]);
            ADASRateInterpolator *SCD_interper = adas->GetSCD(Zs[iz]);
            
            real_t dWi = 0;
            real_t Li = 0;
            real_t Bi = 0;
            
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                len_t indZ = ionHandler->GetIndex(iz,Z0);
                if(Zs[iz]==1 && opacity_modes[iz]==OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){
		            for (len_t i = 0; i < NCells; i++){
		                Li =  amjuel->getIonizLossLyOpaque(Z0, n_cold[i], T_cold[i]);
	                    Li += amjuel->getRecRadLyOpaque(Z0, n_cold[i], T_cold[i]);
				        
				        Bi = 0;
				        // Binding energy rate term
				        if(Z0>0){       // Recombination gain
				        	dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0-1);
				            Bi -= dWi * amjuel->getRecLyOpaque(Z0, n_cold[i], T_cold[i]);
        	                if(includePRB)
	                            Li+=bremsPrefactor*sqrt(T_cold[i])*Z0*Z0*(1 + bremsRel1*T_cold[i]/Constants::mc2inEV);
			            }

                        if (Z0 == 0)
                            diffWeights[NCells*indZ + i] = neutralScale * (Li+Bi);
                        else
                            diffWeights[NCells*indZ + i] = Li+Bi;
	                }
	            }else{
		            for (len_t i = 0; i < NCells; i++){
		                Li =  PLT_interper->Eval(Z0, n_cold[i], T_cold[i]);
		                if (includePRB)
		                    Li += PRB_interper->Eval(Z0, n_cold[i], T_cold[i]);
		                Bi = 0;
		                if(Z0>0 && includePRB)
		                    Bi -= dWi * ACD_interper->Eval(Z0, n_cold[i], T_cold[i]);
		                if(Z0<Zs[iz]){
		                    dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0);
		                    Bi += dWi * SCD_interper->Eval(Z0, n_cold[i], T_cold[i]);
		                }

                        if (Z0 == 0)
                            diffWeights[NCells*indZ + i] = neutralScale * (Li+Bi);
                        else
                            diffWeights[NCells*indZ + i] = Li+Bi;
		            }
                }
            }
        }
        if(!includePRB){ //bremsstrahlung contribution
            for(len_t i=0; i<NCells; i++){
                for(len_t iz = 0; iz<nZ; iz++)
                    for(len_t Z0=0; Z0<=Zs[iz]; Z0++){
                        len_t indZ = ionHandler->GetIndex(iz,Z0);
                        real_t dIonTerm = Z0*Z0;
                        real_t dRelativisticCorrection = bremsRel1*dIonTerm;

                        if (Z0 == 0)
                            diffWeights[NCells*indZ + i] += neutralScale * bremsPrefactor*sqrt(T_cold[i])*(dIonTerm + dRelativisticCorrection*T_cold[i]/Constants::mc2inEV);
                        else
                            diffWeights[NCells*indZ + i] += bremsPrefactor*sqrt(T_cold[i])*(dIonTerm + dRelativisticCorrection*T_cold[i]/Constants::mc2inEV);
                    }
            }
        }
    } else if(derivId == id_ncold){
        for(len_t iz = 0; iz<nZ; iz++){
            ADASRateInterpolator *PLT_interper = adas->GetPLT(Zs[iz]);
            ADASRateInterpolator *PRB_interper = adas->GetPRB(Zs[iz]);
            ADASRateInterpolator *ACD_interper = adas->GetACD(Zs[iz]);
            ADASRateInterpolator *SCD_interper = adas->GetSCD(Zs[iz]);
            
            real_t dWi = 0;
            real_t dLi = 0;
            real_t dBi = 0;
            
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                len_t indZ = ionHandler->GetIndex(iz,Z0);
                if(Zs[iz]==1 && opacity_modes[iz]==OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){
		            for (len_t i = 0; i < NCells; i++){
		                dLi =  amjuel->getIonizLossLyOpaque_deriv_n(Z0, n_cold[i], T_cold[i]);
	                    dLi += amjuel->getRecRadLyOpaque_deriv_n(Z0, n_cold[i], T_cold[i]);
				        
				        dBi = 0;
				        // Binding energy rate term
				        if(Z0>0){       // Recombination gain
				        	dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0-1);
				            dBi -= dWi * amjuel->getRecLyOpaque_deriv_n(Z0, n_cold[i], T_cold[i]);
				            if(includePRB)
				                dLi+=0.5*bremsPrefactor/sqrt(T_cold[i])*Z0*Z0*(1 + 3.0*bremsRel1*T_cold[i]/Constants::mc2inEV);
			            }
                        
                        if (Z0 == 0)
                            diffWeights[i] += neutralScale * n_i[indZ*NCells + i]*(dLi+dBi);
                        else
                            diffWeights[i] += n_i[indZ*NCells + i]*(dLi+dBi);
	                }
                }else{                
		            for (len_t i = 0; i < NCells; i++){
		                real_t dLi = PLT_interper->Eval_deriv_n(Z0, n_cold[i], T_cold[i]);
		                if (includePRB)
		                    dLi += PRB_interper->Eval_deriv_n(Z0, n_cold[i], T_cold[i]);
		                real_t dBi = 0;
		                if(Z0>0 && includePRB)
		                    dBi -= dWi * ACD_interper->Eval_deriv_n(Z0, n_cold[i], T_cold[i]);
		                if(Z0<Zs[iz]){
		                    dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0);
		                    dBi += dWi * SCD_interper->Eval_deriv_n(Z0, n_cold[i], T_cold[i]);
		                }

                        if (Z0 == 0)
                            diffWeights[i] += neutralScale * n_i[indZ*NCells + i]*(dLi+dBi);
                        else
                            diffWeights[i] += n_i[indZ*NCells + i]*(dLi+dBi);
		            }
                }
            }
        }
        if(!includePRB) //bremsstrahlung contribution
            for(len_t i=0; i<NCells; i++)
                diffWeights[i] += bremsPrefactor*sqrt(T_cold[i])*bremsRel2*T_cold[i]/Constants::mc2inEV;
    } else if (derivId == id_Tcold){
        for(len_t iz = 0; iz<nZ; iz++){
            ADASRateInterpolator *PLT_interper = adas->GetPLT(Zs[iz]);
            ADASRateInterpolator *PRB_interper = adas->GetPRB(Zs[iz]);
            ADASRateInterpolator *ACD_interper = adas->GetACD(Zs[iz]);
            ADASRateInterpolator *SCD_interper = adas->GetSCD(Zs[iz]);
            
            real_t dWi = 0;
            real_t dLi = 0;
            real_t dBi = 0;
            
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                len_t indZ = ionHandler->GetIndex(iz,Z0);
                if(Zs[iz]==1 && opacity_modes[iz]==OptionConstants::OPACITY_MODE_GROUND_STATE_OPAQUE){
		            for (len_t i = 0; i < NCells; i++){
		                dLi =  amjuel->getIonizLossLyOpaque_deriv_T(Z0, n_cold[i], T_cold[i]);
		                dLi += amjuel->getRecRadLyOpaque_deriv_T(Z0, n_cold[i], T_cold[i]);
				        
				        dBi = 0;
				        // Binding energy rate term
				        if(Z0>0){       // Recombination gain
				        	dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0-1);
				            dBi -= dWi * amjuel->getRecLyOpaque_deriv_T(Z0, n_cold[i], T_cold[i]);
			            }

                        if (Z0 == 0)
                            diffWeights[i] += neutralScale * n_i[indZ*NCells + i]*(dLi+dBi);
                        else
                            diffWeights[i] += n_i[indZ*NCells + i]*(dLi+dBi);
	                }
                }else{ 
		            for (len_t i = 0; i < NCells; i++){
		                real_t dLi = PLT_interper->Eval_deriv_T(Z0, n_cold[i], T_cold[i]);
		                if (includePRB)
		                    dLi += PRB_interper->Eval_deriv_T(Z0, n_cold[i], T_cold[i]);
		                real_t dBi = 0;
		                if(Z0>0 && includePRB)
		                    dBi -= dWi * ACD_interper->Eval_deriv_T(Z0, n_cold[i], T_cold[i]);
		                if(Z0<Zs[iz]){
		                    dWi = Constants::ec * nist->GetIonizationEnergy(Zs[iz],Z0);
		                    dBi += dWi * SCD_interper->Eval_deriv_T(Z0, n_cold[i], T_cold[i]);
		                }
                        
                        if (Z0 == 0)
                            diffWeights[i] += neutralScale * n_i[indZ*NCells + i]*(dLi+dBi);
                        else
                            diffWeights[i] += n_i[indZ*NCells + i]*(dLi+dBi);
		            }
                }
            }
        }
        if(!includePRB){ //bremsstrahlung contribution
            for(len_t i=0; i<NCells; i++){
                real_t ionTerm = ionHandler->GetZeff(i)*ionHandler->GetFreeElectronDensityFromQuasiNeutrality(i);
                real_t relativisticCorrection = bremsRel1*ionTerm + bremsRel2*n_cold[i]; 
                diffWeights[i] += 0.5*bremsPrefactor/sqrt(T_cold[i])*(ionTerm + 3.0*relativisticCorrection*T_cold[i]/Constants::mc2inEV);
            }
        }
    }
}

