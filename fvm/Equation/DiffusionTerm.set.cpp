/**
 * Implementation of a general 'SetXXX()' routine for the
 * FVM DiffusionTerm.
 */

// void DiffusionTerm::SetXXX(f) {
//
    const len_t nr = grid->GetNr();
    len_t offset = 0;

    const real_t
        *dr   = grid->GetRadialGrid()->GetDr(),
        *dr_f = grid->GetRadialGrid()->GetDr_f();
    
    // Iterate over interior radial grid points
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = grid->GetMomentumGrid(ir);

        const len_t
            np1 = mg->GetNp1(),
            np2 = mg->GetNp2();

        const real_t
            *Vp     = grid->GetVp(ir),
            *Vp_fr  = grid->GetVp_fr(ir),
            *Vp_fr1 = grid->GetVp_fr(ir+1),
            *Vp_f1  = grid->GetVp_f1(ir),
            *Vp_f2  = grid->GetVp_f2(ir),
            *dp1    = mg->GetDp1(),
            *dp2    = mg->GetDp2(),
            *dp1_f  = mg->GetDp1_f(),
            *dp2_f  = mg->GetDp2_f();

        for (len_t j = 0; j < np2; j++) {
            // Do not set terms in the negative trapped region where the 
            // distribution is mirrored and Vp=0
            if(grid->IsNegativePitchTrappedIgnorableCell(ir,j))
                continue; 

            for (len_t i = 0; i < np1; i++) {
                real_t S;

                /////////////////////////
                // RADIUS
                /////////////////////////
                
                #define X(K,V) f((K),i,j,(V))
                
                // Phi^(r)_{k-1/2}
                if (ir > 0) {
                    // XXX: Here, we explicitly assume that the momentum grids are
                    // the same at all radii, so that p at (ir, i, j) = p at (ir+1, i, j)
                    S = Drr(ir, i, j, drr)*Vp_fr[j*np1+i] / (dr[ir]*dr_f[ir-1]*Vp[j*np1+i]);
                    X(ir-1, -S);
                    X(ir,   +S);
                }

                // Phi^(r)_{k+1/2}
                if (ir < nr-1) {
                    // XXX: Here, we explicitly assume that the momentum grids are
                    // the same at all radii, so that p at (ir, i, j) = p at (ir+1, i, j)
                    S = Drr(ir+1, i, j, drr)*Vp_fr1[j*np1+i] / (dr[ir]*dr_f[ir]*Vp[j*np1+i]);
                    X(ir,   +S);
                    X(ir+1, -S);
                }

                #undef X
                
                #define X(I,J,V) f(ir,(I),(J),(V))
                /////////////////////////
                // MOMENTUM 1/1
                /////////////////////////
                // Phi^(1)_{i-1/2,j}
                if (i > 0) {
                    S = D11(ir, i, j, d11)*Vp_f1[j*(np1+1)+i] / (dp1[i]*dp1_f[i-1]*Vp[j*np1+i]);
                    X(i-1, j, -S);
                    X(i, j,   +S);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1) {
                    S = D11(ir, i+1, j, d11)*Vp_f1[j*(np1+1)+i+1] / (dp1[i]*dp1_f[i]*Vp[j*np1+i]);
                    X(i+1, j, -S);
                    X(i,   j, +S);
                }
                
                /////////////////////////
                // MOMENTUM 2/2
                /////////////////////////
                // Phi^(2)_{i-1/2,j}
                if (j > 0) {
                    S = D22(ir, i, j, d22)*Vp_f2[j*np1+i] / (dp2[j]*dp2_f[j-1]*Vp[j*np1+i]);
                    X(i, j,   +S);
                    X(i, j-1, -S);
                }

                // Phi^(2)_{i+1/2,j}
                if (j < np2-1) {
                    S = D22(ir, i, j+1, d22)*Vp_f2[(j+1)*np1+i] / (dp2[j]*dp2_f[j]*Vp[j*np1+i]);
                    X(i, j+1, -S);
                    X(i, j,   +S);
                }
                
                /////////////////////////
                // MOMENTUM 1/2
                /////////////////////////
                // Phi^(1)_{i-1/2,j}
                if (i > 0 && (j > 0 && j < np2-1)) {
                    S = D12(ir, i, j, d12)*Vp_f1[j*(np1+1)+i] / (dp1[i]*(dp2_f[j]+dp2_f[j-1])*Vp[j*np1+i]);
                    X(i,   j+1, +S);
                    X(i-1, j+1, +S);
                    X(i,   j-1, -S);
                    X(i-1, j-1, -S);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1 && (j > 0 && j < np2-1)) {
                    S = D12(ir,i+1,j, d12)*Vp_f1[j*(np1+1)+i+1]/(dp1[i]*(dp2_f[j]+dp2_f[j-1])*Vp[j*np1+i]);
                    X(i+1, j+1, -S);
                    X(i,   j+1, -S);
                    X(i+1, j-1, +S);
                    X(i,   j-1, +S);
                }
                
                /////////////////////////
                // MOMENTUM 2/1
                /////////////////////////
                // Phi^(2)_{i,j-1/2}
                if (j > 0 && (i > 0 && i < np1-1)) {
                    S = D21(ir,i,j,d21)*Vp_f2[j*np1+i] / (dp2[j]*(dp1_f[i]+dp1_f[i-1])*Vp[j*np1+i]);
                    X(i+1, j-1, +S);
                    X(i+1, j,   +S);
                    X(i-1, j-1, -S);
                    X(i-1, j,   -S);
                }

                // Phi^(2)_{i,j+1/2}
                if (j < np2-1 && (i > 0 && i < np1-1)) {
                    S = D21(ir,i,j+1,d21)*Vp_f2[(j+1)*np1+i] / (dp2[j]*(dp1_f[i]+dp1_f[i-1])*Vp[j*np1+i]);
                    X(i+1, j+1, -S);
                    X(i+1, j,   -S);
                    X(i-1, j+1, +S);
                    X(i-1, j,   +S);
                }
                
                
                #undef X
            }
        }

        offset += np1*np2;
    }

// }
