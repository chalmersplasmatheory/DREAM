/**
 * Set elements for the 'PXiExternalCross' boundary condition.
 */

    const len_t nr = this->grid->GetNr();
    len_t loffset = 0, uoffset = 0;

    // Reset matrix offsets...
    const len_t rowOffset = mat->GetRowOffset();
    const len_t colOffset = mat->GetColOffset();
    mat->SetOffset(rowOffset, 0);

    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *lmg = this->lowerGrid->GetMomentumGrid(ir);
        const MomentumGrid *umg = this->upperGrid->GetMomentumGrid(ir);

        const len_t
            lnp = lmg->GetNp1(), lnxi = lmg->GetNp2(),
            unp = umg->GetNp2(), unxi = umg->GetNp2();

        const real_t
            *lxi_f = lmg->GetP2(),
            *uxi_f = umg->GetP2(),
            *ldp   = lmg->GetDp1(),
            *udp   = umg->GetDp1(),
            *ldxi  = lmg->GetDp2(),
            *udxi  = umg->GetDp2();

        const real_t *Ap  = equation->GetAdvectionCoeff1(ir);
        const real_t *Dpp = equation->GetDiffusionCoeff11(ir);
        const real_t *Dpx = equation->GetDiffusionCoeff12(ir);

        const real_t
            *lVp   = this->lowerGrid->GetVp(ir),
            *uVp   = this->upperGrid->GetVp(ir);

        const real_t *delta1 = equation->GetInterpolationCoeff1(ir);

        if (this->type == TYPE_LOWER) {
            len_t idx  = jl*lnp + lnp-1;
            len_t fidx = loffset + idx;

            real_t S = Ap[idx] / (lVp[idx] * ldp[lnp-1]);
            this->__SetFluxComponent(
                
            );
        } else if (this->type == TYPE_UPPER) {
        } else if (this->type == TYPE_DENSITY) {
        } else
            throw EquationException("PXiExternalCross: Unrecognized boundary type: %d\n", this->type);

        /*for (len_t jl = 0, ju = 0; jl < lnxi; jl++) {
            const len_t idx1 = offset1 + j1*np1 +
                (this->type==TYPE_LOWER ? 0.0 : (np1-1));
            const len_t idx2 = offset2 + j2*np2 +
                (this->type==TYPE_LOWER ? (np2-1) : 0.0);

            // Index for coefficients
            const len_t cidx = (this->type==TYPE_LOWER ? idx2 : idx1);
            // Lower side of boundary
            const len_t lidx = (this->type==TYPE_LOWER ? idx2 : idx1);
            // Upper side of boundary
            const len_t uidx = (this->type==TYPE_LOWER ? idx1 : idx2);

            real_t S = 1.0 / (Vp[j1*np1+np1-1]*dp1[np1-1]*dxi1[j1]);

            // When we're considering the upper grid, we subtract
            // the flux from the lower grid.
            if (this->type == TYPE_UPPER)
                S = -S;

            // Patch together the flux across the boundary
            if (j2 == nxi2)
                j2--;

            do {
                real_t dxiBar =
                    min(xi1_f[j1+1], xi2_f[j2+1]) - max(xi1_f[j1], xi2_f[j1]);

                // Advection
                f(idx1, lidx, Ap[cidx]*(1-delta1[ir][cidx]) * S);
                f(idx1, uidx, Ap[cidx]*delta1[ir][cidx] * S);

                if (xi2_f[j2] < xi1_f[j1])
                    j2++;
            } while (j2 < nxi2 && xi2_f[j2] < xi1_f[j1]);
        }*/

        loffset += lnp*lnxi;
        uoffset += unp*unxi;
    }

    // Restore matrix offset
    mat->SetOffset(rowOffset, colOffset);

