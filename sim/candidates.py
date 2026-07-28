"""Initial actuator candidates for the 16 mm-OD envelope.

Shared assumptions:
- Envelope: outer radius <= 8.0 mm, total length incl. 5 mm stroke <= 35 mm.
- Magnet: N52, D8 (r=4) cylinders in commodity lengths — cheap and strong.
- Radial build: magnet r=4.0, guide/bobbin wall + clearance 0.75 -> coil
  r_in=4.75; coil r_out=7.75 (air core, 0.25 housing) or 7.5 (steel shell
  r 7.5..8.0).
- Magnet z_center = 0 is mid-stroke; sweep z in [-2.5, +2.5].
- Coil base NI chosen so ni_scale ~ O(1) at 0.8 N; sign convention:
  push-pull pairs get opposite-sign base NI so scale=+1 drives both usefully.

Run:  python sim/candidates.py [name ...]   (inside `nix develop`)
"""

from __future__ import annotations

import sys

from femm import Coil, Magnet, Model, Steel
from analyze import evaluate, summarize

Z_SWEEP = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
NI_SCALES_LINEAR = [-1.0, -0.5, 0.0, 0.5, 1.0]
NI_SCALES_STEEL = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]

CANDIDATES: dict[str, Model] = {}


def _add(m: Model) -> Model:
    CANDIDATES[m.name] = m
    return m


# A: dual-coil push-pull, air core. Coils flank the magnet's mid-plane;
# opposite polarity so the axial-gradient forces on both poles add.
_add(Model(
    name="aircore_dual",
    magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
    coils=[
        Coil("lo", 4.75, 7.75, -12.5, -0.5, ni=-300.0),
        Coil("hi", 4.75, 7.75, 0.5, 12.5, ni=+300.0),
    ],
))

# B: single coil, air core — cheapest possible build. Magnet rides with its
# top pole inside the coil, bottom pole below the coil mouth.
_add(Model(
    name="aircore_single",
    magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
    coils=[
        Coil("main", 4.75, 7.75, -2.0, 14.0, ni=300.0),
    ],
))

# C: dual coil + steel flux-return shell (0.5 mm wall). Expect higher N/W,
# nonzero cogging, F-vs-I nonlinearity; needs eddy-lag AC check before trusting
# the 10 ms slew spec.
_add(Model(
    name="steel_shell_dual",
    magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
    coils=[
        # 0.05 mm air gap to the shell: avoids coincident collinear segments
        Coil("lo", 4.75, 7.45, -12.5, -0.5, ni=-300.0),
        Coil("hi", 4.75, 7.45, 0.5, 12.5, ni=+300.0),
    ],
    steels=[Steel("shell", 7.5, 8.0, -13.0, 13.0)],
))

# D: long magnet (D8x20), two short coils each centered on a pole at
# mid-stroke — each pole stays coupled to its own coil across the stroke.
_add(Model(
    name="longmag_dual",
    magnet=Magnet(radius=4.0, length=20.0, z_center=0.0),
    coils=[
        Coil("lo", 4.75, 7.75, -13.5, -6.5, ni=-300.0),
        Coil("hi", 4.75, 7.75, 6.5, 13.5, ni=+300.0),
    ],
))


def main(names: list[str]):
    from femm import run_sweep
    for name in names or list(CANDIDATES):
        model = CANDIDATES[name]
        scales = NI_SCALES_STEEL if model.steels else NI_SCALES_LINEAR
        rows = run_sweep(model, Z_SWEEP, scales, tag="v0")
        print(summarize(evaluate(model, rows)))
        print()


if __name__ == "__main__":
    main(sys.argv[1:])

# NOTE (2026-07-28): exploration is done. The three surviving finalist designs
# (opp24c3, w_o3_m12_r, st_w03) live in sim/designs.py as importable builders
# with their verified metadata; results/EXPLORATION.md has the full digest.
# The candidates above are the round-0 starting points, kept for history.
