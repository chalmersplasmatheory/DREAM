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
            // Do not set terms in the negative trapped region where the 
            // distribution is mirrored and Vp=0
            if(grid->IsNegativePitchTrappedIgnorableCell(ir,j))
                continue; 
            bool isNegativeTrappedRadial = grid->IsNegativePitchTrappedIgnorableRadialFluxCell(ir,j);

            for (len_t i = 0; i < np1; i++) {
                real_t 
                    S_i, // advection coefficient on left-hand face of the cell
                    S_o; // advection coefficient on right-hand face of the cell
                real_t *delta;
                const len_t ij = j*np1 + i;
                const real_t vp_ij = Vp[ij];
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

                // Trapping BC: even if the cell is not ignorable, it may still 
                // be such that the radial flux should be mirrored 
                const real_t fr_i = Fr(ir, i, j, fr);
                const real_t fr_o = Fr(ir+1, i, j, fr);
                if( !isNegativeTrappedRadial && ( fr_i || fr_o ) ) {
                    S_i = fr_i *  Vp_fr[ij] / (vp_ij * dr[ir]);
                    S_o = fr_o * Vp_fr1[ij] / (vp_ij * dr[ir]);

                    if(set==JACOBIAN_SET_LOWER){
                        S_i *= 1 - deltaRadialFlux[ir];
                        S_o = 0;
                    } else if(set==JACOBIAN_SET_CENTER) {
                        S_i *= deltaRadialFlux[ir];
                        S_o *= 1 - deltaRadialFlux[ir];
                    } else if(set==JACOBIAN_SET_UPPER) {
                        S_i = 0;
                        S_o *= deltaRadialFlux[ir];
                    }

                    delta = deltar->GetCoefficient(ir,i,j,interp_mode);
                    len_t n;
                    len_t kmin = deltar->GetKmin(ir, &n);
                    len_t kmax = deltar->GetKmax(ir,nr);
                    // Phi^(r)_{ir-1/2,i,j}: Flow into the cell from the "left" r face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(k, -S_i * delta[n]);

                    delta = deltar->GetCoefficient(ir+1,i,j,interp_mode);
                    kmin = deltar->GetKmin(ir+1, &n);   
                    kmax = deltar->GetKmax(ir+1, nr);
                    // Phi^(r)_{ir+1/2,i,j}: Flow out from the cell to the "right" r face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(k,  S_o * delta[n]);
                }
                #undef X
                
                if(set==JACOBIAN_SET_LOWER || set==JACOBIAN_SET_UPPER)
                    continue;

                /////////////////////////
                // MOMENTUM 1
                /////////////////////////
                #define X(I,J,V) f(ir,(I),(J),(V))

                const real_t f1_i = F1(ir, i, j, f1);
                const real_t f1_o = F1(ir, i + 1, j, f1);
                const real_t f1_zero = F1PSqAtZero(ir,j,f1pSqAtZero);
                if(
                    (mg->GetP1_f(i)==0 && f1_zero) || f1_i || f1_o
                ) {
                    real_t VpDp1 = vp_ij*dp1[i];
                    if(mg->GetP1_f(i)==0){
                        // treats singular p=0 point separately
                        const real_t *VpOverP2AtZero = grid->GetVpOverP2AtZero(ir);
                        S_i = f1_zero * VpOverP2AtZero[j] / VpDp1;
                    } else 
                        S_i = f1_i * Vp_f1[ij + j] / VpDp1;
                    S_o = f1_o * Vp_f1[ij + j + 1] / VpDp1;

                    delta = delta1->GetCoefficient(ir,i,j,interp_mode);
                    len_t n;
                    len_t kmin = delta1->GetKmin(i, &n);
                    len_t kmax = delta1->GetKmax(i, np1);
                    // Phi^(1)_{ir,i-1/2,j}: Flow into the cell from the "left" p1 face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(k, j, -S_i * delta[n]);

                    delta = delta1->GetCoefficient(ir,i+1,j,interp_mode);
                    kmin = delta1->GetKmin(i+1, &n);
                    kmax = delta1->GetKmax(i+1, np1);
                    // Phi^(1)_{ir,i+1/2,j}: Flow out from the cell to the "right" p1 face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(k, j,  S_o * delta[n]);
                }
                /////////////////////////
                // MOMENTUM 2
                /////////////////////////

                const real_t f2_i = F2(ir, i, j, f2);
                const real_t f2_o = F2(ir, i, j + 1, f2);
                if(f2_i || f2_o){
                    real_t VpDp2 = vp_ij*dp2[j];
                    S_i = f2_i * Vp_f2[ij] / VpDp2;
                    S_o = f2_o * Vp_f2[ij + np1] / VpDp2;

                    delta = delta2->GetCoefficient(ir,i,j,interp_mode);
                    len_t n;
                    len_t kmin = delta2->GetKmin(j, &n);
                    len_t kmax = delta2->GetKmax(j, np2);
                    // Phi^(2)_{ir,i,j-1/2}: Flow into the cell from the "left" p2 face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(i, k, -S_i * delta[n]);

                    delta = delta2->GetCoefficient(ir,i,j+1,interp_mode);
                    kmin = delta2->GetKmin(j+1, &n);
                    kmax = delta2->GetKmax(j+1, np2);
                    // Phi^(2)_{ir,i,j+1/2}: Flow out from the cell to the "right" p2 face
                    for(len_t k = kmin; k <= kmax; k++, n++)
                        X(i, k,  S_o * delta[n]);
                }
                #undef X
            }
        }

        offset += np1*np2;
    }
