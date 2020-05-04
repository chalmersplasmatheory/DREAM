/**
 *
 * Common implementation of the 'SetMatrixElements()' and 'SetVectorElements()'
 * methods of the 'IonRateEquation' class.
 */

    const len_t Nr = this->grid->GetNr();
    const len_t Z  = this->ions->GetZ(iIon);
    const real_t *n_cold = this->unknowns->GetUnknownData(id_n_cold);
    const len_t ionidx = this->ions->GetIndex(iIon, Z0);

    for (len_t ir = 0; ir < Nr; ir++) {
        // (I_i^(j-1) n_cold + Imp_i^(j-1)) * n_i^(j-1)
        if (Z0 > 0) {
            /*mat->SetElement(
                rOffset+ir, rOffset-Nr+ir,
                Ion[ionidx-1][ir] * n_cold[ir] + Imp[ionidx-1][ir]
            );*/
            NI(-1, Ion[ionidx-1][ir] * n_cold[ir] + Imp[ionidx-1][ir]);
        }

        // -(I_i^(j) n_cold + Imp_i^(j)) * n_i^(j)
        /*mat->SetElement(
            rOffset+ir, rOffset+ir,
            -(Ion[ionidx][ir] * n_cold[ir] + Imp[ionidx][ir])
        );*/
        NI(0, -(Ion[ionidx][ir] * n_cold[ir] + Imp[ionidx][ir]));

        // R_i^(j+1) n_cold * n_i^(j+1)
        if (Z0 < Z) {
            /*mat->SetElement(
                rOffset+ir, rOffset+Nr+ir,
                Rec[ionidx+1][ir] * n_cold[ir]
            );*/
            NI(+1, Rec[ionidx+1][ir] * n_cold[ir]);
        }

        // -R_i^(j) n_cold * n_i^(j-1)
        /*mat->SetElement(
            rOffset+ir, rOffset+ir,
            -Rec[ionidx][ir] * n_cold[ir]
        );*/
        NI(0, -Rec[ionidx][ir] * n_cold[ir]);
    }
