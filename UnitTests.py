# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 19:09:05 2023

@author: Alex
"""
import os
import Martinoid

print("Martinoid.__version__:", Martinoid.__version__)

os.makedirs("UnitTests", exist_ok=True)
os.chdir("UnitTests")

Sequence = "Na-Na-Na"
peptoid = Martinoid.Martinoid(sequence = Sequence, 
                    N_ter_charged=True, 
                    C_ter_charged=False, 
                    N_ter_protect=False,
                    C_ter_protect=False,
                    Linear=False, 
                    Helical=True)

