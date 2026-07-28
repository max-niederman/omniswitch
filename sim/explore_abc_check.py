"""Open-boundary (ABC radius) sensitivity check for aircore_dual.

Runs z=[-2.5, 0], ni_scales=[1.0] at abc_radius 40 / 60 (baseline) / 120 and
reports % force change vs radius. Unique model name per radius so results/
files don't collide; tag='abc'.
"""

import dataclasses
import sys

sys.path.insert(0, "/home/max/Projects/hw/omniswitch/sim")

from femm import run_sweep
from candidates import CANDIDATES

base = CANDIDATES["aircore_dual"]
Z = [-2.5, 0.0]
SCALES = [1.0]

results = {}
for r in (60.0, 40.0, 120.0):
    m = dataclasses.replace(base, name=f"aircore_dual_r{int(r)}", abc_radius=r)
    rows = run_sweep(m, Z, SCALES, tag="abc")
    results[r] = {row["z"]: row["Fz"] for row in rows}
    print(f"abc_radius={r}: " + ", ".join(
        f"z={z:+.1f} Fz={fz:.6f} N" for z, fz in sorted(results[r].items())))

print("\n% change vs abc_radius=120 (reference) and vs 60 (baseline):")
for r in (40.0, 60.0, 120.0):
    for z in Z:
        fz = results[r][z]
        d120 = 100.0 * (fz - results[120.0][z]) / abs(results[120.0][z])
        d60 = 100.0 * (fz - results[60.0][z]) / abs(results[60.0][z])
        print(f"R={r:5.0f} z={z:+.1f}: Fz={fz:.6f} N  "
              f"dv120={d120:+.4f}%  dv60={d60:+.4f}%")
