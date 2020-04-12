/**
 * Implementation of the "set elements" method for the
 * 'MomentQuantity' class.
 */

//void SetElements(X)

    const len_t nr  = fGrid->GetNr();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        // TODO Should probably be something slightly different
        // (we probably don't want to include the 'r' part?)
        const real_t *Vp = fGrid->GetVp(ir);
        const FVM::MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
        const real_t
            *dp1 = mg->GetDp1(),
            *dp2 = mg->GetDp2();

        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();
        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                len_t idx = j*np1 + i;
                X(ir, i, j, integrand[offset+idx] * Vp[idx] * dp1[i] * dp2[j]);
            }
        }

        offset += np1*np2;
    }
