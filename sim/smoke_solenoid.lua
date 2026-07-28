-- Physics smoke test: thin-shell air-core solenoid, axisymmetric.
-- L = 100 mm, mean winding radius R = 5 mm (shell 4.75..5.25 mm), NI = 1000 A-t.
--
-- Analytic checks (thin finite solenoid, Bz uniform in bore at center):
--   Bz_center = mu0*(NI/L) * (L/2)/sqrt((L/2)^2 + R^2) = 12.504 mT
--   flux through r=4mm disk at z=0: Bz_center * pi * (4mm)^2 = 6.285e-7 Wb
--   self-force on winding = 0 by symmetry
--
-- IMPORTANT: do NOT use point evaluation (mo_getb / mo_getpointvalues) — it
-- hangs FEMM under Wine. Contour/block integrals and circuit properties work.

OUT = "Z:/home/max/Projects/hw/omniswitch/results/smoke_solenoid.txt"
FEM = "Z:/home/max/Projects/hw/omniswitch/results/smoke_solenoid.fem"

newdocument(0)                                -- 0 = magnetics
mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)

r1 = 4.75
r2 = 5.25
zh = 50

mi_addnode(r1, -zh)
mi_addnode(r2, -zh)
mi_addnode(r2, zh)
mi_addnode(r1, zh)
mi_addsegment(r1, -zh, r2, -zh)
mi_addsegment(r2, -zh, r2, zh)
mi_addsegment(r2, zh, r1, zh)
mi_addsegment(r1, zh, r1, -zh)

mi_getmaterial("Air")
mi_getmaterial("Copper")
mi_addcircprop("coil", 1000, 1)               -- 1000 A x 1 turn = 1000 A-turns

mi_addblocklabel(5, 0)                        -- winding, group 1
mi_selectlabel(5, 0)
mi_setblockprop("Copper", 1, 0, "coil", 0, 1, 1)
mi_clearselected()

mi_makeABC(7, 250, 0, 0, 0)                   -- open boundary
mi_addblocklabel(2, 0)                        -- air (single connected region)
mi_selectlabel(2, 0)
mi_setblockprop("Air", 1, 0, "", 0, 0, 0)
mi_clearselected()

mi_saveas(FEM)
mi_analyze()
mi_loadsolution()

mo_addcontour(0, 0)
mo_addcontour(4, 0)
flux = mo_lineintegral(0)                     -- flux through disk r=0..4mm
mo_clearcontour()

mo_groupselectblock(1)
vol = mo_blockintegral(10)
fz = mo_blockintegral(19)
mo_clearblock()

analytic = 6.285e-7
h = openfile(OUT, "w")
write(h, "flux_disk4mm_Wb,", flux, "\n")
write(h, "flux_analytic_Wb,", analytic, "\n")
write(h, "rel_err,", (flux - analytic) / analytic, "\n")
write(h, "winding_vol_m3,", vol, "\n")
write(h, "Fz_selfforce_N,", fz, "\n")
closefile(h)

mo_close()
quit()
