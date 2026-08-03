# Piano-feel force requirements (2026-08-02 literature study)

What force envelope a switch needs to replicate a grand-piano key, from the
measurement literature (3-agent sweep; primary sources read: Kinoshita et al.
JASA 121:2959 (2007), Goebl/Bresin/Galembo JASA 118:1154 (2005), Parlitz
J.Biomech 31:1063 (1998), Hirschkorn Waterloo thesis 2004, Somma Acta Acustica
8 (2024), Timmermans et al. Machines 8(4):76 (2020), Igrec *Pianos Inside
Out*, Stanwood; Askenfelt & Jansson 1990–92 via two independent full-text
quotations).

## Static curve (technician regulation values, key front)

- Down weight 0.45–0.55 N (50 gf bass → 46 gf treble), up weight 0.20–0.25 N
  → friction hysteresis ±0.10–0.15 N; balance weight ~0.35 N.
- Hold-down after strike ~0.6 N.
- Key dip **9.5–11.0 mm** (nominal 10).
- Let-off ("escapement") notch: occupies the last 1–2 mm (starts ~6–7.5 mm of
  a 9.5–10 mm dip, resolves by ~8–9 mm). Magnitude quasi-static on a
  regulated, leaded action: **+0.15–0.25 N** bump then a drop; the only
  instrumented slow-moving press (unleaded demonstrator, 0.01 m/s) shows up
  to +1.2 N — speed- and regulation-dependent, flagged uncertain.
- Damper pickup adds ~+0.2–0.3 N from ~3–4 mm.

## Dynamic playing forces (measured at the key)

- Legato pp barely exceeds 0.5 N; struck pp ~2–3 N.
- Staccato peaks: ~8 N (p), ~15 N (mf), up to 50 N (ff); 4/10 expert pianists
  exceeded 60 N (Kinoshita). Legato ≈ 1/3 of staccato at equal dynamic.
- ≥80 % of the impulse lands before key-bottom at soft dynamics (~65 % at ff
  — the rest is keybed reaction, i.e. supplied by a physical end stop).
- Key velocity ≤ ~1 m/s even at forte; descent 25 ms (f) – 160 ms (p);
  3–5 ms impact spike on struck touch.

## Where the big forces come from: reflected inertia

m_eff = m_key' + α²·m_hammer, α (key→hammer velocity ratio) ≈ 5–6, hammer
11–13 g bass → 4–5 g treble ⇒ **equivalent mass at the key ≈ 0.30–0.40 kg
bass/mid, 0.15–0.20 kg treble** (derived from Hirschkorn's measured inertias;
no published total in grams exists — treat as a range). Hammer ≈ 70–85 % of
it. At playing accelerations (tens of m/s²) this reproduces the measured
8–50 N. After let-off the hammer decouples — the reflected inertia *drops*
mid-stroke, a rendered event, not a static curve feature.

Haptic-key precedents: Timmermans/Dehez (UCLouvain) sized their voice coil to
**"50 N in <20 ms"** for full realism; Oboe's MIKEY proves the opposite end —
statics offloaded passively, actuator renders only dynamics. Passive weighted
digital actions "miss escapement, check and repetition, the most significant
haptic effects" (Timmermans) — i.e. the *signature* events are sub-1.5 N.

## Implications for st_w10

| Tier | Needs | st_w10 status |
|---|---|---|
| Static curve, friction, let-off, damper, hold | ≤ ~1.5 N, 0.15–0.3 N features | ✓ (0.8 N cont., ~2.6 N transient stock; 3.8 N w/ 1.5 A driver) |
| pp–mp inertia rendering (0.15–0.4 kg, ≤10 m/s²) | 1.5–4 N transient | ✓ marginal with 1.5 A driver |
| mf–ff staccato realism | 15–50 N, 5–20 ms | ✗ 4–13× over; s≈7–13 drive — outside verified linearity (s=2) and demag envelope; needs a scaled-up actuator (piano key pitch is 23.5 mm — more OD room than 16) |
| Key dip | 9.5–11 mm | ✗ 5 mm stroke. A travel-amplifying lever (~2:1) restores dip but **halves** force at the finger (0.8→0.4 N cont.) — water-cooled continuous (~2–3.5 N) absorbs this; peak becomes ~1.9 N |

Verdict: convincing piano at soft/moderate dynamics incl. every signature
discrete event is in envelope with the 1.5 A driver; concert-grade ff inertia
is an order of magnitude out and was never this design's target. Firmware
must render inertia actively (double-differentiated Hall position → F=m_eff·a
feedforward) — the 4.7 g mover supplies ~1–3 % of it physically.
