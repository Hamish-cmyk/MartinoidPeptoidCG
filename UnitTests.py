# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 19:09:05 2023

@author: Alex
"""
from pathlib import Path
import os
import Martinoid

print("Martinoid.__version__:", Martinoid.__version__)

os.makedirs("UnitTests", exist_ok=True)
os.chdir("UnitTests")
this_directory = Path(".")

Sequence = "Na-Na-Na"
peptoid = Martinoid.Martinoid(sequence = Sequence, 
                    N_ter_charged=True, 
                    C_ter_charged=False, 
                    N_ter_protect=False,
                    C_ter_protect=False,
                    Linear=False, 
                    Helical=True)



itp = (this_directory / f"{Sequence}.itp").read_text()

assert "1	Qd	1	Nx	BB	 1  1.0 ; " in itp
assert "3	P3	3	Nx	BB	 3  0.0 ; " in itp
assert "1	2	3	2	100.0	700.0 ; BB-BB-BB" in itp


print("PASSED")


Sequence = "Na-Na-Na"
peptoid = Martinoid.Martinoid(sequence = Sequence, 
                    N_ter_charged=True, 
                    C_ter_charged=False, 
                    N_ter_protect=False,
                    C_ter_protect=False,
                    Linear=False, 
                    Helical=False)



itp = (this_directory / f"{Sequence}.itp").read_text()

assert "1	Qd	1	Nx	BB	 1  1.0 ; " in itp
assert "3	P3	3	Nx	BB	 3  0.0 ; " in itp
assert "1	2	3	2	138.0	50.0 ; BB-BB-BB" in itp


print("PASSED")