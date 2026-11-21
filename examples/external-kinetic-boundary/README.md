# External kinetic boundary under field reversal

Run from the DREAM repository root with a built `dreami` executable and a
Python environment with DREAM's dependencies and Matplotlib.

```sh
python examples/external-kinetic-boundary/reversal.py reversal-demo \
    --binary build/iface/dreami
```

The script runs a cylindrical runaway-grid case twice. A passing electron
packet starts near the outer momentum face, with runaway sources and radial
transport disabled. The field drives the packet toward the face, then reverses
from +1 to −1 V/m. Both runs use the same initial packet. The control uses
DREAM's `BC_F_0` boundary and the two-way run enables
`setExternalKineticBoundary()`.

The script checks the change in each pitch inventory against its recorded
boundary flux. It also checks the total density budget and the signed
contribution to the runaway current.

| Output | What it shows |
| --- | --- |
| `density-budget.png` and `density-budget.pdf` | Prescribed field, resolved and external densities, and accumulated outward and return transfer |
| `distribution-contours.png` and `distribution-contours.pdf` | A positive-pitch view near the outer face at the start, at peak external density, and after reversal, on a shared color scale |
| `summary.json` | Peak external density, transferred fractions, final resolved density, and normalized balance errors |
| `two-way/` and `legacy/` | Settings, native output, and logs for both runs |

Add `--analyze-only` to regenerate the figures from the saved runs.
