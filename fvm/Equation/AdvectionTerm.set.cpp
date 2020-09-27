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
            for (len_t i = 0; i < np1; i++) {
                real_t 
                    S_i, // advection coefficient on left-hand face of the cell
                    S_o; // advection coefficient on right-hand face of the cell
                real_t *delta;
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


                S_i = Fr(ir,   i, j, fr) *  Vp_fr[j*np1+i] / (Vp[j*np1+i] * dr[ir]);
                S_o = Fr(ir+1, i, j, fr) * Vp_fr1[j*np1+i] / (Vp[j*np1+i] * dr[ir]);

                delta = deltar->GetCoefficient(ir,i,j,interp_mode);
                // Phi^(r)_{ir-1/2,i,j}: Flow into the cell from the "left" r face
                for(len_t n, k = deltar->GetKmin(ir, &n); k <= deltar->GetKmax(ir,nr); k++, n++)
                    X(k, -S_i * delta[n]);

                delta = deltar->GetCoefficient(ir+1,i,j,interp_mode);
                // Phi^(r)_{ir+1/2,i,j}: Flow out from the cell to the "right" r face
                for(len_t n, k = deltar->GetKmin(ir+1, &n); k <= deltar->GetKmax(ir+1,nr); k++, n++)
                    X(k,  S_o * delta[n]);
                
                #undef X
                
                /////////////////////////
                // MOMENTUM 1
                /////////////////////////
                #define X(I,J,V) f(ir,(I),(J),(V))

                if(mg->GetP1_f(i)==0){
                    // treats singular p=0 point separately
                    const real_t *VpOverP2AtZero = grid->GetVpOverP2AtZero(ir);
                    S_i = F1PSqAtZero(ir,j,f1pSqAtZero) * VpOverP2AtZero[j] / (Vp[j*np1]*dp1[i]);
                } else 
                    S_i = F1(ir, i, j, f1) * Vp_f1[j*(np1+1) + i] / (Vp[j*np1+i]*dp1[i]);
                S_o = F1(ir, i+1, j, f1) * Vp_f1[j*(np1+1) + i+1] / (Vp[j*np1+i]*dp1[i]);

                delta = delta1->GetCoefficient(ir,i,j,interp_mode);
                // Phi^(1)_{ir,i-1/2,j}: Flow into the cell from the "left" p1 face
                for(len_t n, k = delta1->GetKmin(i,&n); k <= delta1->GetKmax(i,np1); k++, n++)
                    X(k, j, -S_i * delta[n]);
                // Phi^(1)_{ir,i+1/2,j}: Flow out from the cell to the "right" p1 face
                delta = delta1->GetCoefficient(ir,i+1,j,interp_mode);
                for(len_t n, k = delta1->GetKmin(i+1,&n); k <= delta1->GetKmax(i+1,np1); k++, n++)
                    X(k, j,  S_o * delta[n]);
                

                /////////////////////////
                // MOMENTUM 2
                /////////////////////////


                S_i = F2(ir, i, j,   f2) * Vp_f2[j*np1+i]     / (Vp[j*np1+i]*dp2[j]);
                S_o = F2(ir, i, j+1, f2) * Vp_f2[(j+1)*np1+i] / (Vp[j*np1+i]*dp2[j]);

                delta = delta2->GetCoefficient(ir,i,j,interp_mode);
                // Phi^(2)_{ir,i,j-1/2}: Flow into the cell from the "left" p2 face
                for(len_t n, k = delta2->GetKmin(j,&n); k <= delta2->GetKmax(j,np2); k++, n++)
                    X(i, k, -S_i * delta[n]);

                delta = delta2->GetCoefficient(ir,i,j+1,interp_mode);
                // Phi^(2)_{ir,i,j+1/2}: Flow out from the cell to the "right" p2 face
                for(len_t n, k = delta2->GetKmin(j+1,&n); k <= delta2->GetKmax(j+1,np2); k++, n++)
                    X(i, k,  S_o * delta[n]);
                
                #undef X
            }
        }

        offset += np1*np2;
    }

// }
