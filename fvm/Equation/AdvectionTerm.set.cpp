/**
 * Implementation of the general 'SetXXX' routine. It is
 * made specific with the help of the 'f' macro.
 */

// void AdvectionTerm::SetXXX(f) {

    const len_t nr = grid->GetNr();
    len_t offset = 0;

    const real_t
        *dr = grid->GetRadialGrid()->GetDr();

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
            *dp2    = mg->GetDp2();

        for (len_t j = 0; j < np2; j++) {
            // Evaluate flux in first point
            //real_t S = F1[j*(np1+1) + 1] * h2_f1[j*(np1+1) + 1] * h3_f1[j*(np1+1) + 1] / dp1[1];

            for (len_t i = 0; i < np1; i++) {
                /////////////////////////
                // RADIUS
                /////////////////////////
                #define X(K,V) f((K),i,j,V)
                // XXX: Here we assume that the momentum grid is the same at all
                // radial points.
                //
                // OTHERWISE...:
                // This term is pretty difficult, since we need to evaluate the
                // coefficients in different radial grid points, which in general
                // have different momentum grids and thus require interpolation
                // across grids. For this application, we should be able to assume
                // that momentum grids at all radii use the same coordinate systems.
                // In general, it is much more difficult, though (so it would require
                // a bit more thinking if we wanted to interpolate generally between
                // two different momentum grids)

                // Phi^(r)_{ir-1/2,i,j}
                if (ir > 0) {
                    real_t S = Fr(ir, i, j, fr) * Vp_fr[j*np1+i] / (Vp[j*np1+i] * dr[ir]);
                    X(ir-1, -S * (1-deltar[ir][j*np1 + i]));
                    X(ir,   -S * deltar[ir][j*np1 + i]);
                }

                // Phi^(r)_{ir+1/2,i,j}
                if (ir < nr-1) {
                    real_t S = Fr(ir+1, i, j, fr) * Vp_fr1[j*np1+i] / (Vp[j*np1+i] * dr[ir]);
                    X(ir,   S * (1-deltar[ir+1][j*np1 + i]));
                    X(ir+1, S * deltar[ir+1][j*np1 + i]);
                }
                
                #undef X
                
                /////////////////////////
                // MOMENTUM 1
                /////////////////////////
                #define X(I,J,V) f(ir,(I),(J),(V))

                // Phi^(1)_{i-1/2,j}
                if (i > 0) {
                    real_t S = F1(ir, i, j, f1) * Vp_f1[j*(np1+1) + i] / (Vp[j*np1+i]*dp1[i]);
                    X(i-1, j,-S * (1-delta1[ir][j*np1 + i]));
                    X(i,   j,-S * delta1[ir][j*np1 + i]);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1) {
                    real_t S = F1(ir, i+1, j, f1) * Vp_f1[j*(np1+1) + i+1] / (Vp[j*np1+i]*dp1[i]);
                    X(i,   j, S * (1-delta1[ir][j*np1 + i+1]));
                    X(i+1, j, S * delta1[ir][j*np1 + i+1]);
                }

                /////////////////////////
                // MOMENTUM 2
                /////////////////////////
                // Phi^(2)_{i,j-1/2}
                if (j > 0) {
                    real_t S = F2(ir, i, j, f2) * Vp_f2[j*np1+i] / (Vp[j*np1+i]*dp2[j]);
                    X(i, j-1,-S * (1-delta2[ir][j*np1+i]));
                    X(i, j,  -S * delta2[ir][j*np1+i]);
                }

                // Phi^(2)_{i,j+1/2}
                if (j < np2-1) {
                    real_t S = F2(ir, i, j+1, f2) * Vp_f2[(j+1)*np1+i] / (Vp[j*np1+i]*dp2[j]);
                    X(i, j,   S * (1-delta2[ir][(j+1)*np1+i]));
                    X(i, j+1, S * delta2[ir][(j+1)*np1+i]);
                }

                #undef X
            }
        }

        offset += np1*np2;
    }

// }
