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
        const real_t *xi_f = mg->GetP2_f();

        // Determine xi integration limits
        // (note that since xi=0 may lie anywhere within a grid cell,
        // we must properly adjust the first/last cell length to
        // accurately correspond to the actual positive/negative part
        // of the cell)
        len_t j0 = 0, jmax = np2;
        real_t dxi_first = dp2[0], dxi_last = dp2[np2-1];
        switch (xiMode) {
            case XI_MODE_NEG: {
                for (len_t j = 1; j < np2+1; j++) {
                    if (xi_f[j]*xi_f[j-1] < 0) {
                        jmax = j;
                        dxi_last = -xi_f[j-1];
                        break;
                    }
                }
            } break;
            case XI_MODE_POS: {
                for (len_t j = 1; j < np2+1; j++) {
                    if (xi_f[j]*xi_f[j-1] < 0) {
                        j0 = j-1;
                        dxi_first = xi_f[j];
                        break;
                    }
                }
            } break;

            default: break;
        }

        for (len_t i = 0; i < np1; i++){
            real_t envelope = ThresholdEnvelope(ir,i);
            real_t common = envelope * dp1[i] / VpVol;

            X(j0*np1+i, common * Vp[j0*np1+i] * dxi_first);
            for (len_t j = j0+1; j < jmax-1; j++){
                len_t idx = j*np1 + i;
                X(idx, common * Vp[idx] * dp2[j])
            }
            X((jmax-1)*np1+i, common * Vp[(jmax-1)*np1+i] * dxi_last);
        }
        ApplyX(ir)
        offset += np1*np2;
    }
