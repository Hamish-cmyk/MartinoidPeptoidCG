# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:57:27 2023

@author: rkb19187
"""

import numpy as np
import matplotlib.pyplot as plt

bond = np.genfromtxt("angle5_a0.xvg")

print(bond)

plt.plot(bond[:,0], bond[:,1], 'g', lw=2)
plt.ylim(-2,50)
plt.xlim(70,)
plt.xlabel("Angle ($^{\circ}$)")
plt.ylabel("Energy")

fig = plt.gcf()
fig.savefig("images/CustomAngle.png", dpi=300)