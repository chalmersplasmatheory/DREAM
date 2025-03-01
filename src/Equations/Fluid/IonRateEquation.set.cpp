/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    
    //const len_t ionidx = this->ions->GetIndex(iIon, Z0);
    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // I_i^(j-1) n_cold * n_i^(j-1)
            if (Z0 > 0)
                NI(-1, Ion[Z0-1][ir] * n_cold[ir], posIoniz);
		    
            // -I_i^(j) n_cold * n_i^(j)
            NI(0, -Ion[Z0][ir] * n_cold[ir], negIoniz);
        }
        // R_i^(j+1) n_cold * n_i^(j+1)
        if (Z0 < Z)
            NI(+1, Rec[Z0+1][ir] * n_cold[ir], posRec);

        // -R_i^(j) n_cold * n_i^(j-1)
        NI(0, -Rec[Z0][ir] * n_cold[ir], negRec);
    }

	// Charge-exchange
	if (cxIons.size() > 0) {
		const len_t nIons = this->ions->GetNZ();

		const real_t *Wi = unknowns->GetUnknownData(id_W_i);
		const real_t *Ni = unknowns->GetUnknownData(id_N_i);

		for (len_t cxIon : cxIons) {
			ADASRateInterpolator *ccdIon = GetCCD(iIon);

			// Iterate over radii
			for (len_t ir = 0; ir < Nr; ir++) {
				real_t Ti;
				if (Ni[iIon*Nr+ir] <= 0) Ti = 0;
				else Ti = 2.0/3.0 * Wi[iIon*Nr+ir] / (DREAM::Constants::ec*Ni[iIon*Nr+ir]);

				if (iIon == cxIon) {
					real_t ni = ions->GetIonDensity(ir, iIon, 1);
					real_t Rcx_i = ccdIon->Eval(0, ni, Ti);

					for (len_t kIon = 0; kIon < nIons; kIon++) {
						// CX with itself does not change particle content
						if (kIon == iIon)
							continue;

						len_t Zk = ions->GetZ(kIon);
						const len_t ionOffset = ions->GetIndex(kIon, 0);

						ADASRateInterpolator *ccd = GetCCD(kIon);

						// Iterate over k charge states
						for (len_t kZ0 = 1; kZ0 <= Zk; kZ0++) {
							real_t nk = ions->GetIonDensity(ir, kIon, kZ0);
							len_t idx = kIon*Nr + ir;

							real_t Tk;
							if (Ni[idx] <= 0) Tk = 0;
							else Tk = 2.0/3.0 * Wi[idx] / (DREAM::Constants::ec*Ni[idx]);

							real_t Rcx = ccd->Eval(kZ0-1, nk, kZ0);

							if (Z0 == 0) {
								// Apply to neutral hydrogen
								NI(0, -Rcx * nions[(ionOffset+kZ0)*Nr+ir], negCX);

								// Extra term for when both species are hydrogenic
								if (Zk == 1)
									NI(1, +Rcx_i * nions[(ionOffset+0)*Nr+ir], posCX);
							} else /*if (Z0 == 1)*/ {
								// Apply to neutral hydrogen (Z0-1 = 0)
								NI(-1, Rcx * nions[(ionOffset+kZ0)*Nr+ir], posCX);

								// Extra term for when both species are hydrogenic
								if (Zk == 1)
									NI(0, -Rcx_i * nions[(ionOffset+0)*Nr+ir], negCX);
							}
						}
					}
				} else if (Z0 < Z) {
					real_t nZ0 = ions->GetIonDensity(ir, iIon, Z0+1);
					real_t Rcx_i = ccdIon->Eval(Z0, nZ0, Ti);

					for (len_t kIon = 0; kIon < nIons; kIon++) {
						// Only CX which involves at least one hydrogenic
						// species counts
						if (ions->GetZ(kIon) != 1)
							continue;

						const len_t kOffset = ions->GetIndex(kIon, 0);
						NI(+1, Rcx_i * nions[kOffset*Nr+ir], posCX);
					}
				}

				// Negative CX term
				if (Z != 1 && Z0 >= 1) {
					real_t nZ0 = ions->GetIonDensity(ir, iIon, Z0);

					real_t Rcx_i = ccdIon->Eval(Z0-1, nZ0, Ti);
					for (len_t kIon = 0; kIon < nIons; kIon++) {
						if (ions->GetZ(kIon) != 1)
							continue;

						const len_t kOffset = ions->GetIndex(kIon, 0);
						NI(0, -Rcx_i * nions[kOffset*Nr+ir], negCX);
					}
				}
			}
		}
	}
