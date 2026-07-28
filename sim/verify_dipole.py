"""Analytic cross-check: FEMM force on a small magnet vs point-dipole model.

Geometry (synthetic, ignores package envelope on purpose):
  * one thin coil, r 9.5..10.5 mm, z 0..1 mm, NI = 100 A-turns
    -> quasi-current-loop of radius a = 10 mm centered at z = 0.5 mm
  * small N52 magnet r = 1.5 mm, l = 3 mm, swept z_center in [-20, -15, -10]

Analytic model:
  m  = Br * V / mu0                      (A m^2, axial dipole)
  Bz(z) = mu0 * NI * a^2 / (2 (a^2+z^2)^1.5)   (on-axis loop field, z from loop plane)
  dBz/dz = -3 mu0 NI a^2 z / (2 (a^2+z^2)^2.5)
  F  = m * dBz/dz  evaluated at the magnet center (z measured from loop plane).

Expectation: few-% agreement at 15-20 mm separation, monotonically worse closer
(magnet finite size + off-axis field curvature over the magnet volume).
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from femm import MU0, Coil, Magnet, Model, run_sweep  # noqa: E402

A_LOOP = 10.0e-3        # loop radius, m (coil mean radius)
Z_LOOP = 0.5e-3         # loop plane, m (axial center of the thin coil)
NI = 100.0              # ampere-turns

MAG_R = 1.5e-3          # m
MAG_L = 3.0e-3          # m
BR = 1.43               # T (N52)

Z_POSITIONS = [-20.0, -15.0, -10.0]   # magnet z_center, mm


def dipole_moment() -> float:
    vol = math.pi * MAG_R ** 2 * MAG_L
    return BR * vol / MU0            # A m^2


def dbz_dz(z_from_loop: float) -> float:
    """Axial gradient of on-axis loop field, T/m. z_from_loop in m."""
    a2 = A_LOOP ** 2
    return -3.0 * MU0 * NI * a2 * z_from_loop / (2.0 * (a2 + z_from_loop ** 2) ** 2.5)


def analytic_force(z_center_mm: float) -> float:
    z = z_center_mm * 1e-3 - Z_LOOP
    return dipole_moment() * dbz_dz(z)


def main() -> None:
    model = Model(
        name="verify_dipole",
        magnet=Magnet(radius=1.5, length=3.0, z_center=Z_POSITIONS[0],
                      br=BR, mu_r=1.05),
        coils=[Coil("loop", r_in=9.5, r_out=10.5, z_bot=0.0, z_top=1.0,
                    ni=NI)],
        abc_radius=60.0,
        mesh_air=0.25,   # fine mesh in magnet + coil; force is only ~1 mN
    )
    # The interaction force (~mN) is far below the magnet self-force mesh
    # noise in the stress tensor integral.  Since F is linear in I, the
    # antisymmetric difference (F(+NI) - F(-NI))/2 cancels the
    # current-independent self-force error exactly; the ni_scale=0 row
    # reports that error (noise floor) directly.
    rows = run_sweep(model, Z_POSITIONS, [-1.0, 0.0, 1.0], tag="dipole")

    by_z: dict[float, dict[float, float]] = {}
    for row in rows:
        by_z.setdefault(row["z"], {})[row["ni_scale"]] = row["Fz"]

    m = dipole_moment()
    print(f"dipole moment m = {m:.6e} A m^2")
    print(f"{'z_mm':>7} {'F_femm_N':>13} {'F_dipole_N':>13} {'err_%':>8} "
          f"{'noise_N(ni=0)':>14}")
    for z in sorted(by_z, key=lambda v: v):
        f = by_z[z]
        f_femm = 0.5 * (f[1.0] - f[-1.0])
        f_an = analytic_force(z)
        err = 100.0 * (f_femm - f_an) / f_an
        print(f"{z:7.1f} {f_femm:13.6e} {f_an:13.6e} {err:8.2f} "
              f"{f[0.0]:14.3e}")


if __name__ == "__main__":
    main()
