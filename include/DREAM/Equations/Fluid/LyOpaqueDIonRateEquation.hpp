#ifndef _DREAM_EQUATION_LY_OPAQUE_D_ION_RATE_EQUATION_HPP
#define _DREAM_EQUATION_LY_OPAQUE_D_ION_RATE_EQUATION_HPP

#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "DREAM/AMJUEL.hpp"

namespace DREAM {
    class LyOpaqueDIonRateEquation : public IonRateEquation {
	private:
    	AMJUEL *amjuel;
    public:
    	LyOpaqueDIonRateEquation(
    		FVM::Grid* g, IonHandler* ions, const len_t iIon,
    		FVM::UnknownQuantityHandler* unknowns, bool addFluidIonization, bool addFluidJacobian, bool isAbl, AMJUEL* amjuel
    	):IonRateEquation(g,ions,iIon,nullptr,unknowns,addFluidIonization, addFluidJacobian,isAbl), amjuel(amjuel){}
 
    	virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override{
    		const len_t Nr = this->grid->GetNr();

			real_t *T = unknowns->GetUnknownData(id_T_cold);
			real_t *n = unknowns->GetUnknownData(id_n_cold);

			real_t eps = sqrt(std::numeric_limits<real_t>::epsilon());
			// Iterate over charge state (0 ... Z)
			for (len_t i = 0; i < Nr; i++){
				real_t hn = eps*(1 + n[i]);
				real_t hT = eps*(1 + T[i]);
				for (len_t Z0 = 0; Z0 <= Zion; Z0++){
					Rec[Z0][i]         = amjuel->getRecLyOpaque(Z0, n[i], T[i]);
					PartialNRec[Z0][i] = (amjuel->getRecLyOpaque(Z0, n[i]+hn, T[i]) - amjuel->getRecLyOpaque(Z0, n[i], T[i]))/hn;
					PartialTRec[Z0][i] = (amjuel->getRecLyOpaque(Z0, n[i], T[i]+hT) - amjuel->getRecLyOpaque(Z0, n[i], T[i]))/hT;
					Ion[Z0][i]         = 0;
					PartialNIon[Z0][i] = 0;
					PartialTIon[Z0][i] = 0;
				}
			}
			// if not covered by the kinetic ionization model, set fluid ionization rates
			if(addFluidIonization || addFluidJacobian)
				for (len_t i = 0; i < Nr; i++){
					real_t hn = eps*(1 + n[i]);
					real_t hT = eps*(1 + T[i]);
					for (len_t Z0 = 0; Z0 <= Zion; Z0++){
						Ion[Z0][i]         = amjuel->getIonizLyOpaque(Z0, n[i], T[i]);
						PartialNIon[Z0][i] = (amjuel->getIonizLyOpaque(Z0, n[i]+hn, T[i]) - amjuel->getIonizLyOpaque(Z0, n[i], T[i]))/hn;
						PartialTIon[Z0][i] = (amjuel->getIonizLyOpaque(Z0, n[i], T[i]+hT) - amjuel->getIonizLyOpaque(Z0, n[i], T[i]))/hT;
					}
				}
		}
	};
}

#endif/*_DREAM_EQUATION_LY_OPAQUE_D_ION_RATE_EQUATION_HPP*/
