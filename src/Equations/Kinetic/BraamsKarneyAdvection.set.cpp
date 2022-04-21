/**
 * Implementation of the general 'SetXXX' routine. It is
 * made specific with the help of the 'f' macro.
 */

// void AdvectionTerm::SetXXX(f) {

    const len_t nr = grid->GetNr();
    len_t offset = 0;

    // const real_t
    //     *dr = grid->GetRadialGrid()->GetDr();

    // Iterate over interior radial grid points
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);

        const len_t
            np1 = mg->GetNp1(),
            np2 = mg->GetNp2();

        const real_t
            *Vp     = grid->GetVp(ir),
            // *Vp_fr  = grid->GetVp_fr(ir),
            // *Vp_fr1 = grid->GetVp_fr(ir+1),
            *Vp_f1  = grid->GetVp_f1(ir),
            *Vp_f2  = grid->GetVp_f2(ir),
            *dp1    = mg->GetDp1(),
            *dp2    = mg->GetDp2();

        for (len_t j = 0; j < np2; j++) {
            // Do not set terms in the negative trapped region where the 
            // distribution is mirrored and Vp=0
            if(grid->IsNegativePitchTrappedIgnorableCell(ir,j))
                continue; 
            //bool isNegativeTrappedRadial = grid->IsNegativePitchTrappedIgnorableRadialFluxCell(ir,j);

            for (len_t i = 0; i < np1; i++) {
                real_t 
                    S_i, // advection coefficient on left-hand face of the cell
                    S_o; // advection coefficient on right-hand face of the cell

                real_t *delta;

                /////////////////////////
                // MOMENTUM 1
                /////////////////////////
                #define X(I,J,V) f(ir,(I),(J),(V))

                if(
                    (mg->GetP1_f(i)==0 && F1PSqAtZero(ir,j,f1pSqAtZero)) ||
                    F1(ir, i, j, dfs[k2][k1].df1) || (k1 > 0 && F1(ir, i+1, j, dfs[k2][k1 - 1].df1))
                ) {
                    real_t VpDp1 = Vp[j*np1+i]*dp1[i];
                    if(mg->GetP1_f(i)==0){
                        // treats singular p=0 point separately
                        const real_t *VpOverP2AtZero = grid->GetVpOverP2AtZero(ir);
                        S_i = F1PSqAtZero(ir,j,f1pSqAtZero) * VpOverP2AtZero[j] / VpDp1;
                    } else 
                        S_i = F1(ir, i, j, dfs[k2][k1].df1) * Vp_f1[j*(np1+1) + i] / VpDp1;
					if (k1 > 0)
						S_o = F1(ir, i+1, j, dfs[k2][k1 - 1].df1) * Vp_f1[j*(np1+1) + i+1] / VpDp1;
					else
						S_o = 0;

					delta = delta1->GetCoefficient(ir,i,j,interp_mode);
                    // Phi^(1)_{ir,i-1/2,j}: Flow into the cell from the "left" p1 face
                    for(len_t n, k = delta1->GetKmin(i,&n); k <= delta1->GetKmax(i,np1); k++, n++)
                        X(k, j, -S_i * delta[n]);
                    // Phi^(1)_{ir,i+1/2,j}: Flow out from the cell to the "right" p1 face
                    delta = delta1->GetCoefficient(ir,i+1,j,interp_mode);
                    for(len_t n, k = delta1->GetKmin(i+1,&n); k <= delta1->GetKmax(i+1,np1); k++, n++)
                        X(k, j,  S_o * delta[n]);
                }
                /////////////////////////
                // MOMENTUM 2
                /////////////////////////

                if(F2(ir, i, j, dfs[k2][k1].df2) || (k2 > 0 && F2(ir, i, j+1, dfs[k2 - 1][k1].df2))){
                    real_t VpDp2 = Vp[j*np1+i]*dp2[j];
                    S_i = F2(ir, i, j, dfs[k2][k1].df2) * Vp_f2[j*np1+i]     / VpDp2;
					if (k2 > 0)
						S_o = F2(ir, i, j+1, dfs[k2 - 1][k1].df2) * Vp_f2[(j+1)*np1+i] / VpDp2;
					else
						S_o = 0;

					delta = delta2->GetCoefficient(ir,i,j,interp_mode);
                    // Phi^(2)_{ir,i,j-1/2}: Flow into the cell from the "left" p2 face
                    for(len_t n, k = delta2->GetKmin(j,&n); k <= delta2->GetKmax(j,np2); k++, n++)
                        X(i, k, -S_i * delta[n]);

                    delta = delta2->GetCoefficient(ir,i,j+1,interp_mode);
                    // Phi^(2)_{ir,i,j+1/2}: Flow out from the cell to the "right" p2 face
                    for(len_t n, k = delta2->GetKmin(j+1,&n); k <= delta2->GetKmax(j+1,np2); k++, n++)
                        X(i, k,  S_o * delta[n]);
                }
                #undef X
            }
        }

        offset += np1*np2;
    }
