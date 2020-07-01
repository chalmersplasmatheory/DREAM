/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    //const len_t ionidx = this->ions->GetIndex(iIon, Z0);

    for (len_t ir = 0; ir < Nr; ir++) {
        // (I_i^(j-1) n_cold + Imp_i^(j-1)) * n_i^(j-1)
        if (Z0 > 0) {
            NI(-1, Ion[Z0-1][ir] );
        }

        // -(I_i^(j) n_cold + Imp_i^(j)) * n_i^(j)
        NI(0, -Ion[Z0][ir]);

        // R_i^(j+1) n_cold * n_i^(j+1)
        if (Z0 < Z) {
            NI(+1, Rec[Z0+1][ir] );
        }

        // -R_i^(j) n_cold * n_i^(j-1)
        NI(0, -Rec[Z0][ir]);
    }
