# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:39:12 2023

@author: rkb19187
"""

import Martinoid

def readin(fname):
    f = open(fname, 'r', errors="ignore")
    content = f.read()
    return content

peptoid = Martinoid.Martinoid(sequence = "Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv", N_ter_charged=True, SS="Helical")
assert "90.0	400.0" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "1	Qd" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "20	Nda	8	Nv	BB" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "100.0	700.0 ; BB-BB-BB" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")


peptoid = Martinoid.Martinoid(sequence = "Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv", N_ter_charged=False,  SS="linear")
assert "0.0	20.0	1 ; BB-BB-BB-BB " in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "1	Nda	1	Nf	BB	 1  0.0 ;" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "20	Nda	8	Nv	BB" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "138.0	50.0 ; BB-BB-BB" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")

peptoid = Martinoid.Martinoid(sequence = "Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv", C_ter_charged=True, C_ter_protect=False, SS="")
assert "0.0	20.0	1 ; BB-BB-BB-BB " not in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "90.0	400.0" not in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "20	Qa	8	Nv	BB	 20  -1.0 ; " in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "1	Qd	1	Nf	BB	 1  1.0 ; " in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")

peptoid = Martinoid.Martinoid(sequence = "Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv", N_ter_charged=False, N_ter_protect=True, SS="helical")
assert "1	Na	1	Nf	BB	 1  0.0 ;" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "90.0	400.0" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
assert "20	Nda	8	Nv	BB	 20  0.0 ;" in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")



# =============================================================================
# # This one should FAIL
# peptoid = Martinoid.Martinoid(sequence = "Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv", N_ter_charged=False, C_ter_charged=True, SS="")
# assert "0.0	20.0	1 ; BB-BB-BB-BB " not in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
# assert "90.0	400.0" not in readin("Nf-Nfe-Nq-Nm-NmO-Nv-Nv-Nv.itp")
# =============================================================================


