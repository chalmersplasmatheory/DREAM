# Two-way external kinetic boundary

Run from the DREAM repository root with a built `dreami` and a
Python environment with DREAM's dependencies and Matplotlib.

```sh
python examples/external-kinetic-boundary/reversal.py two-pole-demo \
    --binary build/iface/dreami
```

A dilute, isotropic seed starts near zero momentum. A prescribed electric
field first accelerates it to positive parallel momentum, reverses to drive
it through zero to negative parallel momentum, then reverses again and
brings it back near zero. Each excursion reaches just beyond one unit of
parallel momentum in the collisionless limit. The kinetic grid uses the
nonnegative magnitude p. The two crossings of p=1 are near the positive
and negative ends of the pitch range, respectively.

Both runs use the same initial packet, time steps, and cell widths. Runaway
sources and radial transport are disabled.

| Run | Outer momentum face | Outside p=1 |
| --- | --- | --- |
| `wide/` | p=2, `BC_F_0` | The packet remains on the kinetic grid to within the reported numerical loss through p=2 |
| `cut/` | p=1, `setExternalKineticBoundary()` | The two external densities store particles by pitch sign and limit their return flux |

The script checks the particle balance and the flux into and out of each
external density. It checks that the wide-grid run retains its seed, and
that both packets pass through zero parallel momentum between the excursions
and return near zero by the end of the pulse. The wider grid resolves the
momentum distribution throughout the excursion. The cut grid reconstructs
the return flux from its last accepted boundary cells and tracks only density
outside its kinetic domain. The two runs can therefore differ in their
distribution shapes and return timing. First-order upwind advection also
broadens the packets.

| Output | What it shows |
| --- | --- |
| `density-budget.png` and `density-budget.pdf` | Field pulse, resolved and total densities, and both cut-grid inventories compared with the wide-grid density in the same momentum and pitch ranges |
| `distribution-contours.png` and `distribution-contours.pdf` | Shared logarithmic contours in parallel and perpendicular momentum at five times, with the p=1 face marked |
| `summary.json` | Peak pole populations, final resolved fractions and mean parallel momenta, normalized balance errors, and the maximum wide-grid density loss |
| `cut/` and `wide/` | Settings, DREAM outputs, and native run logs |

Add `--analyze-only` to regenerate the figures from the saved runs.
