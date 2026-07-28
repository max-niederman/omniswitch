# omniswitch actuator — frozen design record (2026-07-28)

Selected design: **st_w10** — single moving magnet, dual coil, 1.0 mm steel
flux-return shell doubling as housing and magnetic shield.

All numbers below are FEMM-verified (axisymmetric magnetostatics + harmonic
eddy study; pipeline validated to <1% against analytic models — see
docs/EXPLORATION.md and the sweep CSVs referenced in sim/converge_shell.py).

## Geometry (mm, radial from axis)

| Part | r | z | Notes |
|---|---|---|---|
| Magnet | 0 – 4.0 | ±6.0 (12 long), travels ±2.5 | D8×12 cylinder, axially magnetized, **N45SH** (N38UH if any 120 °C exposure) |
| Guide gap + bobbin | 4.0 – 4.75 | | PTFE/POM guide surface |
| Coils ×2 | 4.75 – 6.95 | lo −12.5..−0.5, hi 0.5..12.5 | opposed polarity, series, one continuous wire |
| Air | 6.95 – 7.0 | | |
| Shell = housing | 7.0 – 8.0 | ±13 | 1018 steel tube, 16 mm OD |

Package: **28 × 16 mm** incl. 5 mm stroke and end caps (specs: ≤35, ≤16).

## Performance vs requirements

| Requirement | Result | Status |
|---|---|---|
| OD ≤ 16 mm | 16.0 | PASS (at limit) |
| ≥0.8 N everywhere on 5 mm stroke, both directions | gain 2.54–2.82 N per unit s (N52-sim); F ∝ I linear to <0.1%, no saturation to 2× drive | PASS, ~3× headroom |
| Length ≤ 35 mm | 28.0 | PASS |
| ≤5 W continuous at 0.8 N worst point | 1.53 W (20 °C) – 1.95 W (100 °C magnet), fill 0.5, N45SH, incl. cogging | PASS (thermal reality caps ~2–2.5 W continuous; we sit under it) |
| 0→0.8 N < 10 ms at ≤24 V | τ=137 µs, rise ≈49 µs + 5.7 µs shell eddy lag (measured, harmonic solve) | PASS ×~180 |
| <$5 @ 100 units | est. $3.5–4.5 (quotes pending) | PASS (verify) |
| Moving PM, precise force control | F ∝ I exact; cogging 0.161 N clean centering spring (odd in z) → firmware table | PASS |

## Why this one (decision trail)

- Air-core opposed-pair designs (opp24c3, w_o3_m12_r) were 2–4× more
  power-efficient BUT: (a) **crosstalk-fatal at 19.05 mm pitch** — the like-pole
  radial flux band injects 1.1–1.2 N on the neighbor's mover (quadrupole
  cancellation is a far-field story; 19 mm is deep near-field); (b) demag-stressed
  at the facing poles (opp24c3 unusable above ~70 °C even on N38UH); (c) opp24c3
  needs 31 N of permanent retention against pair repulsion (the washer variant
  cuts that to 2.4 N but keeps the crosstalk problem).
- Pot-core stator end washers: REJECTED — 0.4–0.9 N parasitic spring.
- 0.3 mm shell (st_w03): wall saturates at 1.85 T → only ~2× shielding;
  keyboard-matrix analysis says non-viable even with feedforward (residual up
  to 17% of 0.8 N).
- **st_w10's 1.0 mm wall stays at 1.21 T (unsaturated, µ_diff ≈ 325)**:
  measured source-leakage attenuation ×18 (worst neighbor force 25 mN,
  upper bound with unshielded victim), incoming transverse shielding S_in ≈ 23.
  Cost: +26% hold power vs thin shell — absorbed by margin.
- Demag: thicker return path slightly *raises* the magnet operating point.
  N52 disqualified in every design (margin ≤0.90). st_w10 margins (1.5×
  corner-peaking factor, bus-stall reverse overdrive s=1.27):
  N45SH 1.62/1.40/1.17 @ 80/100/120 °C; N38UH 2.20/1.94/1.65.
  Full record: results/demag_margins.csv (regenerate:
  `nix develop -c python sim/qualify_demag.py analyze`).

## Winding + drive (preliminary — re-quote after bobbin design)

~480 turns/coil AWG 33 heavy-build (fill 0.5), two coils series-opposed, one
continuous wire. R ≈ 27 Ω hot; hold at worst case ≈ 270 mA / 7.3 V (duty 0.30
on 24 V); bus-stall ≈ 0.9 A → stall force ≈ 2.6 N. Driver: DRV8871 for
bring-up, AT8870-class clone (~$0.3) for production; current sense via 0.3 Ω
shunt + INA181, mid-ON-time sampling; PWM 25 kHz phase-staggered in groups.
Force calibration: K(z) per-unit factory map on N45SH; temperature-compensate
K by −0.12 %/K of magnet temp using winding resistance as the thermometer.

## Keyboard integration (from the matrix feasibility study)

- 19.05 mm pitch is viable. Feedforward compensation of the neighbor force map
  is **mandatory** (firmware knows all key positions + currents); with st_w10
  leakage the residual is well under the 2–5%-of-0.8 N curve-shaping bar.
- Use **row-stripe magnet polarity** (alternate rows flipped): −40% adversarial
  / −81% deterministic neighbor map vs uniform. Checkerboard is worse — do not use.
- Sensor: neighbor field at own Hall ≈ 0.18 mT (3 µm apparent error) — negligible.
- Power: 104-key all-hold 0.8 N ≈ 62 W worst (24 V / 5 A supply with firmware
  global current budget + e-fuse; realistic typing ~5 W). ≥2 mF bulk +
  4.7 µF/key; 27–30 V TVS for release regen bursts.

## Shell slit: evaluated and REJECTED (2026-07-28)

Question: axial slit(s) in the shell to block eddy currents. Answer: **keep the
plain closed tube** — the slit removes almost nothing here, and costs real
performance. Studies: `sim/slit_shielding.py` (planar 2D transverse-shielding),
`sim/`+`results/st_w10_slitac.csv`, `st_w10_slitside.csv`.

- **The slit doesn't kill the eddy screening in this design.** Induced E is
  azimuthal, but because the coils are driven series-opposed the AC flux is odd
  in z and the net-winding circumferential mode carries ~zero current
  (FEMM: |I_net| ≤ 16 mA while distributed eddies dissipate 6.8 W at solve
  drive). The actual screening is zero-net-winding local loops (azimuthal out
  at one z, back at another, closing axially) — these survive a slit by closing
  along the slit edges. FEMM proof: a shell constrained to zero net current is
  digit-identical to the unslit shell at 20 Hz–50 kHz; a common-mode control
  drive (which this design never produces) shows the slit working (−81% lag).
  Analytic plate-mode estimate: 1 slit keeps ~65% of the lag, 2 slits ~37%.
- Benefit if slit anyway: 5.7 µs → ~2–4 µs lag (spec margin already 1750×) and
  1.4 mW → 0.3–0.7 mW PWM ripple loss (vs ~1.4 W hold). Noise.
- **Single-slit cost**: 27–54 mN constant-direction lateral preload on the mover
  (Maxwell stress of the missing wedge; kerf 0.2/0.4 mm) → 4–8 mN
  velocity-sign-dependent rail friction, 40–80% of the 10 mN feel threshold and
  NOT firmware-compensable (unlike cogging). Plus worst-orientation transverse
  shielding halves (S_eff 23 → 11–15).
- **Two symmetric slits**: side-load cancels to a 15–17 mN tolerance residual
  (±0.05 mm kerf, ±2°), BUT both circumferential shunt paths break —
  worst-orientation S_eff collapses to 4.0–5.6 (neighbor force 4.5–6.2 mN,
  approaching feel threshold before feedforward), the slit geometry is
  gap-reluctance-limited (high µ no longer helps), and the housing becomes two
  loose half-shells (~170× torsional derating).
- A slit aligned with the transverse field costs zero shielding — but a keyboard
  victim has orthogonal neighbors on both axes, so the worst orientation binds;
  45°-to-grid slit azimuth recovers ~25–40%.
- If a slit is ever *forced* (manufacturing): ONE slit, ≤0.2–0.3 mm kerf,
  bridged ends (one-piece tube), laser/EDM-cut DOM tube (never roll-formed
  open-seam as the OD-setting housing), azimuth 45° off-grid or toward board
  edge. If eddy screening ever genuinely matters, the fix is a ferrite return,
  not a slit.

## Open items before committing hardware

1. Vendor datasheet + quotes: D8×12 N45SH (margins near threshold use catalog
   Hk=0.9·Hci from memory — sim/magnets.py flags this), steel tube, winding.
2. Bobbin design → exact coil r_in; re-run winding numbers if r_in > 4.75.
3. Prototype validation: force-vs-current-vs-position rig; cogging table
   measurement; two-unit crosstalk measurement at 19.05 mm to confirm the ×18.
4. Position sensor selection + integration (not yet simulated — Hall in bore
   at z ≈ −9.5 assumed by the crosstalk study; shell interaction TBD).
5. End caps: steel caps would add axial shielding + a return path but also
   cogging — currently assumed polymer; simulate if sensor needs quieting.
