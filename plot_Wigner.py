# -*- coding: utf-8 -*-
"""
Example to plot nice Wigner functions

Created on Fri May 17 11:05:01 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import matplotlib.pyplot as plt

N = 20
sys = cv.System(N, 1)

r = 1
sys.apply_SMS(r, 0)
statei = sys.state
#statei = qt.basis(N, 1)   

# Photon subtraction options
t = .9

# Scissors options
# Best  k=0.5 m_aux=0.6
k = .6
m_aux = .2


r_aux = np.arcsinh(np.sqrt(m_aux))
print("Parameter:", r_aux)
p_sc = sys.apply_scissors(k, r_aux)
print("p_SC:", p_sc)
#sys.apply_scissors_exact(k)




sys2 = cv.System(N, 1)

sys2.apply_SMS(r, 0)
p_ps = sys2.apply_photon_subtraction(t, 0)
print("p_PS:", p_ps)
#sys2.apply_scissors_exact_options(k, 0, 'a')


x = np.linspace(-4, 4, 100)
y = np.linspace(-4, 4, 100)

w1 = qt.wigner(sys.state, x, y)
cmap1 = qt.wigner_cmap(w1)

w2 = qt.wigner(statei, x, y)
cmap2 = qt.wigner_cmap(w2)

w3 = qt.wigner(sys2.state, x, y)
cmap3 = qt.wigner_cmap(w3)


fig, axes = plt.subplots(1, 3, figsize=(12,3))
im1 = axes[1].imshow(w1)
axes[1].set_title("Out")

norm = im1.norm

im0 = axes[0].imshow(w2, norm=norm)
axes[0].set_title("Scissors")

im2 = axes[2].imshow(w3, norm=norm)
axes[2].set_title("Photon Subtraction")


#fig.colorbar(im1)
fig.show()




fig, axes = plt.subplots(1, 3, figsize=(12,3))
qt.plot_fock_distribution(statei, fig=fig, ax=axes[0], title="In", unit_y_range=False);
qt.plot_fock_distribution(sys.state, fig=fig, ax=axes[1], title="Scissors", unit_y_range=False);
qt.plot_fock_distribution(sys2.state, fig=fig, ax=axes[2], title="Photon Subtraction", unit_y_range=False);
fig.tight_layout()
plt.show()

