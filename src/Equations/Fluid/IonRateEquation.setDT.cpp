/**
 * Sets temperature jacobian of IonRateEquation
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);

    for (len_t ir = 0; ir < Nr; ir++) {
        if(setIonization){
            // (I_i^(j-1) n_cold
            if (Z0 > 0) 
                NI(-1, PartialTIon[Z0-1][ir] * n_cold[ir]);

            // -(I_i^(j) n_cold + Imp_i^(j)) * n_i^(j)
           NI(0, - PartialTIon[Z0][ir] * n_cold[ir] );
        }
        // R_i^(j+1) n_cold * n_i^(j+1)
        if (Z0 < Z) {
            NI(+1, PartialTRec[Z0+1][ir] * n_cold[ir]);
        }

        // -R_i^(j) n_cold * n_i^(j-1)
        NI(0, -PartialTRec[Z0][ir] * n_cold[ir]);
    }
