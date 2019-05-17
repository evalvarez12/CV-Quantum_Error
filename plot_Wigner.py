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

r = .8
sys.apply_SMD(r, 0)
statei = sys.state
#statei = qt.basis(N, 1)   

k = .1
m_aux = .01
r_aux = np.arcsinh(np.sqrt(m_aux))
print(r_aux)
sys.apply_scissors(k, r_aux)
#sys.apply_scissors_exact(k)




sys2 = cv.System(N, 1)

sys2.apply_SMD(r, 0)
sys2.apply_scissors_inverted(k, r_aux)





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
axes[0].set_title("In")

im0 = axes[2].imshow(w3, norm=norm)
axes[2].set_title("Compare")

fig.colorbar(im1)
fig.show()
