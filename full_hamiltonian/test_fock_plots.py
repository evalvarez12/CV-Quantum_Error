# -*- coding: utf-8 -*-
"""
Plotting some Fock distribtions to compare how calculations over
a finite Hilbert space behave

Created on Mon Apr  1 14:48:35 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import matplotlib.pyplot as plt
import numpy as np


r = 3
print("Mean photon number:", np.sinh(r)**2)

# N = 10 
N = 10
a = qt.basis(N, 0)
Sa = qt.squeeze(N, r)

# N = 20 
N = 20
b = qt.basis(N, 0)
Sb = qt.squeeze(N, r)

# N = 30 
N = 50
c = qt.basis(N, 0)
Sc = qt.squeeze(N, r)


fig, axes = plt.subplots(1, 3, figsize=(12,3))
qt.plot_fock_distribution(Sa*a, fig=fig, ax=axes[0], title="A", unit_y_range=False);
qt.plot_fock_distribution(Sb*b, fig=fig, ax=axes[1], title="B", unit_y_range=False);
qt.plot_fock_distribution(Sc*c, fig=fig, ax=axes[2], title="C", unit_y_range=False);
fig.tight_layout()
plt.show()