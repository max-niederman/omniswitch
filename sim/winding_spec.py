"""Concrete winding + drive-electronics spec for the two finalists
(opp24c3, w_o3_m12_r) on the 24 V bus, honest fill = 0.5, hot copper.

No FEMM runs: works from known geometry + the fine-mesh sweep CSVs
(noise-cancelled force gains, single-turn flux linkages).

Method = the aircore_dual power-model audit method:
  * per-section NI ratios fixed (300/600/300 and 300/115/115), sections in
    series with ONE continuous wire -> same current, turns split by NI ratio;
  * the section with max NI/A_window is binding; fill 0.5 there sets N for a
    given gauge; V_hold = I*R is gauge-set (indep. of N for fixed gauge), so
    the gauge is chosen for hold duty ~0.3 of 24 V at worst stroke position;
  * L_eff (incl. mutuals, correct series signs) from the combined-drive flux
    derivative: L = (dLambda/ds)/(dI/ds), Lambda = sum sgn_k * N_k * lam_k.

Run: nix develop /home/max/Projects/hw/omniswitch -c python sim/winding_spec.py
Writes results/winding_spec.csv
"""

import csv
import math
import os

ROOT = "/home/max/Projects/hw/omniswitch"
RES = os.path.join(ROOT, "results")

V_BUS = 24.0
F_TARGET = 0.8
FILL = 0.5
RHO20 = 1.72e-8
ALPHA_CU = 0.00393          # /K
T_HOT = 75.0                # spec temperature for R_tot
F_PWM = 25e3
RHO_MAG = 7.5               # g/cm^3 (prompted)
RHO_STEEL = 7.87
BR_N52, BR_N45SH, TC_BR = 1.445, 1.32, -0.0012   # sim/magnets.py
BR_SIM = 1.43               # Br actually used in the FEMM models (femm.py)


def rho_cu(T):
    return RHO20 * (1 + ALPHA_CU * (T - 20))


def awg(g):
    d = 0.127e-3 * 92 ** ((36 - g) / 39)   # bare dia, m
    return d, math.pi * d * d / 4


def load(path):
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def gains_and_dlam(rows, coil_names):
    """Noise-cancelled gain (N per unit s) and dlam/ds per coil, per z."""
    byz = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r
    out = {}
    for z, d in byz.items():
        smax, smin = max(d), min(d)
        ds = smax - smin
        gain = (d[smax]["Fz"] - d[smin]["Fz"]) / ds
        dlam = {c: (d[smax][f"flux_{c}"] - d[smin][f"flux_{c}"]) / ds
                for c in coil_names}
        out[z] = (gain, dlam)
    return out


# section: (name, r_in, r_out, z_bot, z_top, ni_base)  sign = winding sense
DESIGNS = {
    "opp24c3": dict(
        csv=os.path.join(RES, "opp24c3f_acd4.csv"),
        sections=[("lo", 4.5, 7.75, -16.5, -8.0, -300.0),
                  ("mid", 4.5, 7.75, -6.5, 6.5, +600.0),
                  ("hi", 4.5, 7.75, 8.0, 16.5, -300.0)],
        ref="mid",
        mover=dict(magnets=2 * [(0.4, 1.2)],       # (r_cm, len_cm) D8x12
                   steel=[], spacer_g=0.02, sleeve_g=0.50),
    ),
    "w_o3_m12_r": dict(
        csv=os.path.join(RES, "w_o3_m12_r_wcgF.csv"),
        sections=[("lo", 4.75, 7.75, -15.5, -5.5, -115.0),
                  ("mid", 4.75, 7.75, -5.0, 5.0, +300.0),
                  ("hi", 4.75, 7.75, 5.5, 15.5, -115.0)],
        ref="mid",
        mover=dict(magnets=2 * [(0.4, 1.2)],
                   steel=[(0.4, 0.2)],             # 2 mm x D8 solid 1018 disc
                   spacer_g=0.0, sleeve_g=0.70),
    ),
}

GAUGE = 34   # chosen below-verified: hold duty ~0.27 (N52) / ~0.30 (N45SH 20C)

out_rows = []


def spec(name, d):
    print(f"\n================ {name} ================")
    secs = d["sections"]
    names = [s[0] for s in secs]
    rows = load(d["csv"])
    g = gains_and_dlam(rows, names)
    zs = sorted(g)
    gain_w = min(g[z][0] for z in zs)
    z_w = [z for z in zs if g[z][0] == gain_w][0]
    gain_c = max(g[z][0] for z in zs)
    s_w = F_TARGET / gain_w

    ref = next(s for s in secs if s[0] == d["ref"])
    ni_ref = abs(ref[5])
    d_bare, a_cu = awg(GAUGE)

    # turns: binding section = ref (max |NI|/A); fill FILL there
    def area(s):
        return (s[2] - s[1]) * (s[4] - s[3]) * 1e-6  # m^2

    A_ref = area(ref)
    N_ref = int(A_ref * FILL / a_cu)
    Ns, senses, fills = {}, {}, {}
    for s in secs:
        Ns[s[0]] = int(round(N_ref * abs(s[5]) / ni_ref))
        senses[s[0]] = "+" if s[5] > 0 else "-"
        fills[s[0]] = Ns[s[0]] * a_cu / area(s)
    N_tot = sum(Ns.values())

    # resistance
    def R(T):
        return sum(rho_cu(T) * 2 * math.pi * (0.5 * (s[1] + s[2]) * 1e-3)
                   * Ns[s[0]] / a_cu for s in secs)

    R20, R75, R100 = R(20), R(T_HOT), R(100)
    wire_len = sum(2 * math.pi * 0.5 * (s[1] + s[2]) * 1e-3 * Ns[s[0]]
                   for s in secs)
    cu_g = wire_len * a_cu * 8.96e6 * 1e-3 * 1000  # m*m^2*kg/m^3 -> g

    # L_eff(z) incl. mutuals: Lambda' = sum sgn*N*dlam/ds ; I' = ni_ref/N_ref
    dIds = ni_ref / Ns[ref[0]]
    Ls = {}
    for z in zs:
        lam_p = sum((1 if s[5] > 0 else -1) * Ns[s[0]] * g[z][1][s[0]]
                    for s in secs)
        Ls[z] = abs(lam_p) / dIds
    z_mid = zs[len(zs) // 2]
    L = Ls[z_mid]
    tau = L / R75

    # hold point (N52, 20C magnet), worst stroke position
    I_hold = s_w * ni_ref / Ns[ref[0]]
    V_hold = I_hold * R75
    duty = V_hold / V_BUS
    P_hold = I_hold * V_hold
    # per-section power split
    P_secs = {s[0]: I_hold ** 2 * rho_cu(T_HOT) * 2 * math.pi
              * (0.5 * (s[1] + s[2]) * 1e-3) * Ns[s[0]] / a_cu for s in secs}

    # force constant K = dF/dI
    K_c = gain_c / dIds
    K_w = gain_w / dIds

    # stall / rise / ripple
    I_stall_hot, I_stall_cold = V_BUS / R75, V_BUS / R20
    F_stall_hot_w = K_w * I_stall_hot
    F_stall_hot_c = K_c * I_stall_hot
    F_stall_cold_c = K_c * I_stall_cold
    t_rise = -tau * math.log(1 - I_hold / I_stall_hot)

    def ripple(f, D):
        return V_BUS * D * (1 - D) / (L * f)

    dI25 = ripple(F_PWM, duty)
    dI25max = ripple(F_PWM, 0.5)
    dI40 = ripple(40e3, duty)

    # N45SH grade change + 80C magnet
    f20 = BR_N45SH / BR_N52
    f80 = BR_N45SH * (1 + TC_BR * (80 - 20)) / BR_N52
    I_hold_45_20 = I_hold / f20
    I_hold_45_80 = I_hold / f80
    V_corner = I_hold_45_80 * R100
    duty_corner = V_corner / V_BUS
    P_corner = I_hold_45_80 * V_corner
    t_rise_corner = -(L / R100) * math.log(1 - I_hold_45_80 / (V_BUS / R100))

    # mover mass + haptics
    mv = d["mover"]
    m_mag = sum(math.pi * r * r * ln * RHO_MAG for r, ln in mv["magnets"])
    m_steel = sum(math.pi * r * r * t * RHO_STEEL for r, t in mv["steel"])
    m = m_mag + m_steel + mv["spacer_g"] + mv["sleeve_g"]
    a3 = 3.0 / (m * 1e-3)
    a08 = 0.8 / (m * 1e-3)
    x08_1ms = 0.5 * a08 * 1e-6 * 1e6      # um in 1 ms
    x3_1ms = 0.5 * a3 * 1e-6 * 1e6

    kv = dict(design=name, gauge_awg=GAUGE, wire_bare_mm=round(d_bare * 1e3, 4),
              wire_area_mm2=round(a_cu * 1e6, 5),
              turns=" / ".join(f"{n}:{senses[n]}{Ns[n]}" for n in names),
              N_total=N_tot, fill_per_sec=" / ".join(
                  f"{n}:{fills[n]:.3f}" for n in names),
              wire_len_m=round(wire_len, 1), copper_g=round(cu_g, 1),
              gain_worst_NperS=round(gain_w, 3), z_worst=z_w,
              gain_center=round(gain_c, 3), s_hold_worst=round(s_w, 4),
              R20=round(R20, 1), R75=round(R75, 1), R100=round(R100, 1),
              L_eff_mH=round(L * 1e3, 2),
              L_range_mH=f"{min(Ls.values())*1e3:.2f}-{max(Ls.values())*1e3:.2f}",
              tau75_us=round(tau * 1e6, 1),
              K_center_NperA=round(K_c, 2), K_worst_NperA=round(K_w, 2),
              I_hold_worst_mA=round(I_hold * 1e3, 1),
              V_hold_worst=round(V_hold, 2), duty_hold=round(duty, 3),
              P_hold_W=round(P_hold, 3),
              P_split_W=" / ".join(f"{n}:{P_secs[n]:.3f}" for n in names),
              I_stall_hot_mA=round(I_stall_hot * 1e3),
              I_stall_cold_mA=round(I_stall_cold * 1e3),
              F_stall_hot_worstpos=round(F_stall_hot_w, 2),
              F_stall_hot_center=round(F_stall_hot_c, 2),
              F_stall_cold_center=round(F_stall_cold_c, 2),
              t_rise_08N_us=round(t_rise * 1e6, 1),
              dI_pp_25k_hold_mA=round(dI25 * 1e3, 1),
              dI_pp_25k_max_mA=round(dI25max * 1e3, 1),
              dI_pp_40k_hold_mA=round(dI40 * 1e3, 1),
              F_ripple_pp_25k_N=round(dI25 * K_c, 3),
              K_scale_N45SH_20C=round(f20, 4),
              K_scale_N45SH_80C=round(f80, 4),
              K_scale_vs_simBr143_20C=round(BR_N45SH / BR_SIM, 4),
              K_scale_vs_simBr143_80C=round(
                  BR_N45SH * (1 + TC_BR * 60) / BR_SIM, 4),
              I_hold_N45SH20_mA=round(I_hold_45_20 * 1e3, 1),
              corner_I_mA=round(I_hold_45_80 * 1e3, 1),
              corner_V=round(V_corner, 2), corner_duty=round(duty_corner, 3),
              corner_P_W=round(P_corner, 3),
              corner_trise_us=round(t_rise_corner * 1e6, 1),
              mover_mass_g=round(m, 2),
              m_magnets_g=round(m_mag, 2), m_steel_g=round(m_steel, 2),
              accel_3N_ms2=round(a3, 0), accel_08N_ms2=round(a08, 1),
              x_1ms_08N_um=round(x08_1ms, 1), x_1ms_3N_um=round(x3_1ms, 0))
    for k, v in kv.items():
        print(f"{k:28s} {v}")
    out_rows.append(kv)
    return kv


for nm, dd in DESIGNS.items():
    spec(nm, dd)

out_csv = os.path.join(RES, "winding_spec.csv")
with open(out_csv, "w", newline="") as f:
    wcsv = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    wcsv.writeheader()
    wcsv.writerows(out_rows)
print(f"\nwrote {out_csv}")

# bifilar-outer option for w_o3_m12_r (outer windows only 19% filled with the
# straight series wind): wind lo+hi 2-in-hand from the same AWG34 spool
# (one joint at the mid boundary) -> outer R halves.
d = DESIGNS["w_o3_m12_r"]
r = out_rows[1]
P_out = sum(float(x.split(":")[1]) for x in r["P_split_W"].split(" / ")
            if not x.startswith("mid"))
P_mid = float([x for x in r["P_split_W"].split(" / ")
               if x.startswith("mid")][0].split(":")[1])
print(f"w_o3_m12_r bifilar-outer option: P_hold {P_mid + P_out/2:.3f} W "
      f"(vs {r['P_hold_W']} W straight series; outer fill 0.38)")
