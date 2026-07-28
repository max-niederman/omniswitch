"""Mesh-convergence check for aircore_dual force pipeline (tag=meshconv)."""

from __future__ import annotations

import dataclasses
import json

from candidates import CANDIDATES
from femm import run_sweep

BASE = CANDIDATES["aircore_dual"]
Z = [0.0, -2.5]

# mesh_air value -> (tag suffix, ni_scales to run)
# noise floor (ni=0) only needed at automesh and finest mesh
RUNS = [
    (0.0, "meshconv_auto", [0.0, 1.0]),
    (0.8, "meshconv_m08", [1.0]),
    (0.4, "meshconv_m04", [1.0]),
    (0.2, "meshconv_m02", [0.0, 1.0]),
]

out = {}
for mesh, tag, scales in RUNS:
    model = dataclasses.replace(BASE, mesh_air=mesh)
    rows = run_sweep(model, Z, scales, tag=tag)
    for r in rows:
        out[(mesh, r["z"], r["ni_scale"])] = r["Fz"]
    print(f"done mesh={mesh} tag={tag}", flush=True)

print(json.dumps([{"mesh": k[0], "z": k[1], "s": k[2], "Fz": v}
                  for k, v in out.items()], indent=1))

# analysis
for z in Z:
    fz = {mesh: out[(mesh, z, 1.0)] for mesh, _, _ in RUNS}
    ref = fz[0.2]
    vals = list(fz.values())
    spread = (max(vals) - min(vals)) / abs(ref) * 100
    print(f"\nz={z}: " + "  ".join(f"mesh={m}: {v:.5f} N" for m, v in fz.items()))
    print(f"  spread (max-min)/|finest| = {spread:.3f} %")
    print(f"  automesh vs finest: {(fz[0.0] - ref) / abs(ref) * 100:+.3f} %")
    for m in (0.0, 0.2):
        if (m, z, 0.0) in out:
            print(f"  noise floor mesh={m}: Fz(ni=0) = {out[(m, z, 0.0)]:.6f} N")
