
#include <cmath>
#include "DREAM/AMJUEL.hpp"

using namespace DREAM;
using namespace std;


/**
 * Lyman-opaque recombination rate coefficient [m^3 s^{-1}]
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecAMJUEL+=recLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return exp(lnRecAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque recombination rate coefficient
 * with respect to the background electron density
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecAMJUEL+=recLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(j>0)
                    derivFactor+=recLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*j*lnnj/lnn/n;
                
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return derivFactor*exp(lnRecAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque recombination rate coefficient
 * with respect to the background electron temperature
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
        	lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecAMJUEL+=recLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(i>0)
                    derivFactor+=recLyOpaque[i*maxDensCoeffsIndex+j]*i*lnTi/lnT/T*lnnj;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return derivFactor*exp(lnRecAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Lyman-opaque recombination radiation coefficient [m^3 W]
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecRadAMJUEL+=recRadLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*exp(lnRecRadAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque recombination radiation coefficient
 * with respect to the background electron density
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecRadAMJUEL+=recRadLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(j>0)
                    derivFactor+=recRadLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*j*lnnj/lnn/n;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque recombination radiation coefficient
 * with respect to the background electron temperature
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getRecRadLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==1){
        real_t lnRecRadAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnRecRadAMJUEL+=recRadLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(i>0)
                    derivFactor+=recRadLyOpaque[i*maxDensCoeffsIndex+j]*i*lnTi/lnT/T*lnnj;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*derivFactor*exp(lnRecRadAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Lyman-opaque ionization rate coefficient [m^3 s^{-1}]
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizAMJUEL+=ionizLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return exp(lnIonizAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }			
}

/**
 * Derivative of the Lyman-opaque ionization rate coefficient
 * with respect to the background electron density
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizAMJUEL+=ionizLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(j>0)
                    derivFactor+=ionizLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*j*lnnj/lnn/n;
                    
                lnnj*=lnn;
            }
            lnTi=lnT;
        }
        return derivFactor*exp(lnIonizAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }			
}

/**
 * Derivative of the Lyman-opaque ionization rate coefficient
 * with respect to the background electron temperature
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizAMJUEL+=ionizLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(i>0)
                    derivFactor+=ionizLyOpaque[i*maxDensCoeffsIndex+j]*i*lnTi/lnT/T*lnnj;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return derivFactor*exp(lnIonizAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }			
}

/**
 * Lyman-opaque ionization loss coefficient, including both the potential energy difference 
 * and line radiation during one effective ionization event (which can include several
 * excitation-deexcitation events before the ion/atom is eventually ionized) [m^3 W]
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*exp(lnIonizLossAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque ionization loss coefficient
 * with respect to the background electron density
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque_deriv_n(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(j>0)
                    derivFactor+=ionizLossLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*j*lnnj/lnn/n;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

/**
 * Derivative of the Lyman-opaque ionization loss coefficient
 * with respect to the background electron temperature
 *
 * Z0: Charge state (either 0 or 1 as these rates only applies to hydrogen isotopes).
 * n:  Background electron density.
 * T:  Background electron temperature.
 */
real_t AMJUEL::getIonizLossLyOpaque_deriv_T(len_t Z0, real_t n, real_t T){
    if(Z0==0){
        real_t lnIonizLossAMJUEL=0;
        real_t derivFactor=0;
        
        real_t lnT=log(T);
        real_t lnn=log(n/1e14);
        real_t lnTi=1.0;
        real_t lnnj=1.0;
        
        for(len_t i=0;i<maxTempCoeffsIndex;i++){
            lnnj=1.0;
            for(len_t j=0;j<maxDensCoeffsIndex;j++){
                lnIonizLossAMJUEL+=ionizLossLyOpaque[i*maxDensCoeffsIndex+j]*lnTi*lnnj;
                if(i>0)
                    derivFactor+=ionizLossLyOpaque[i*maxDensCoeffsIndex+j]*i*lnTi/lnT/T*lnnj;
                    
                lnnj*=lnn;
            }
            lnTi*=lnT;
        }
        return Constants::ec*derivFactor*exp(lnIonizLossAMJUEL)/1e6; // Factor 1e6 converts from cm^3 to m^3
    }else{
        return 0;
    }        			
}

