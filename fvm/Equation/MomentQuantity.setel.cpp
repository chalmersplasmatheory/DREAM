/**
 * Implementation of the "set elements" method for the
 * 'MomentQuantity' class.
 */

    const len_t nr  = fGrid->GetNr();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const real_t *Vp = fGrid->GetVp(ir);
        const real_t VpVol = fGrid->GetVpVol(ir);
        const FVM::MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
        const real_t
            *dp1 = mg->GetDp1(),
            *dp2 = mg->GetDp2();

        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        for (len_t j = 0; j < np2; j++) 
            for (len_t i = 0; i < np1; i++) {
                len_t idx = j*np1 + i;
                X(i, j, ThresholdEnvelope(ir,i,j) * Y(idx) * (Vp[idx] / VpVol) * dp1[i] * dp2[j])
            }
        ApplyX(ir)
        offset += np1*np2;
    }
