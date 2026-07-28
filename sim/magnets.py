"""Magnet grade presets: thermal/demag properties for design selection.

Demagnetization model (second-quadrant, axial magnetization):
  local operating point H_d = (B_d - Br(T)) / (mu0 * mu_r)  [A/m, negative]
  partial demag begins when H_d < -Hk(T)  (intrinsic knee field)
Grade data: typical catalog values (Chinese sintered NdFeB / Sm2Co17).
Tempcos are linear approximations valid ~20..max_temp_C.

Price factors are rough 100-unit relative street prices for a D8x12 cylinder
vs N52 = 1.0 (flag: from memory, verify with quotes before committing BOM).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

MU0 = 4e-7 * math.pi


@dataclass(frozen=True)
class Grade:
    name: str
    br20: float          # T at 20 C
    hci20: float         # intrinsic coercivity A/m at 20 C
    hk20: float          # knee field A/m at 20 C (~0.85-0.95 * Hci for sintered NdFeB)
    tc_br: float         # Br tempco, fraction/K (negative)
    tc_hci: float        # Hci/Hk tempco, fraction/K (negative)
    max_temp_c: float    # catalog max operating temp (Pc ~ 1)
    mu_r: float
    price_factor: float

    def br(self, temp_c: float) -> float:
        return self.br20 * (1 + self.tc_br * (temp_c - 20))

    def hk(self, temp_c: float) -> float:
        return self.hk20 * (1 + self.tc_hci * (temp_c - 20))

    def demag_margin(self, bz_local: float, temp_c: float) -> float:
        """Margin factor = |Hk(T)| / |H_d|. > 1.5 is comfortable, < 1.1 is risky.

        bz_local: local flux density (T, along magnetization) in the magnet
        under worst-case reverse drive, from FEMM sub-block integration.
        """
        h_d = (bz_local - self.br(temp_c)) / (MU0 * self.mu_r)
        if h_d >= 0:
            return float("inf")
        return self.hk(temp_c) / -h_d


GRADES = {
    "N52": Grade("N52", 1.445, 876e3, 0.90 * 876e3, -0.0012, -0.0055, 80, 1.05, 1.0),
    "N45M": Grade("N45M", 1.345, 1114e3, 0.90 * 1114e3, -0.0012, -0.0055, 100, 1.05, 1.1),
    "N45SH": Grade("N45SH", 1.32, 1592e3, 0.90 * 1592e3, -0.0012, -0.0053, 150, 1.05, 1.4),
    "N42SH": Grade("N42SH", 1.30, 1592e3, 0.90 * 1592e3, -0.0012, -0.0053, 150, 1.05, 1.3),
    "N38UH": Grade("N38UH", 1.24, 1990e3, 0.90 * 1990e3, -0.0012, -0.0050, 180, 1.05, 1.5),
    "SmCo26": Grade("SmCo26", 1.05, 1512e3, 0.90 * 1512e3, -0.0003, -0.0020, 300, 1.06, 4.0),
}


def power_scale(from_grade: str, to_grade: str, temp_c: float = 60.0) -> float:
    """Ohmic-power multiplier when swapping grades at equal force (air core):
    P ∝ 1/Br(T)^2."""
    a, b = GRADES[from_grade], GRADES[to_grade]
    return (a.br(temp_c) / b.br(temp_c)) ** 2
