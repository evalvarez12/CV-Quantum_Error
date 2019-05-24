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
sys.apply_SMD(r, 0)
statei = sys.state
#statei = qt.basis(N, 1)   

k = .6
m_aux = .2
r_aux = np.arcsinh(np.sqrt(m_aux))
print(r_aux)
sys.apply_scissors(k, r_aux)
#sys.apply_scissors_exact(k)




sys2 = cv.System(N, 1)

sys2.apply_SMD(r, 0)
sys2.apply_scissors_options(1-k, r_aux, 0, 'a')
#sys2.apply_scissors_exact_options(k, 0, 'a')


sys3 = cv.System(N, 1)

sys3.apply_SMD(r, 0)
sys3.apply_scissors_options(1-k, r_aux, 0, 'b')
#sys3.apply_scissors_exact_options(k, 0, 'b')


sys4 = cv.System(N, 1)

sys4.apply_SMD(r, 0)
sys4.apply_scissors_options(1-k, r_aux, 0, 'c')
#sys4.apply_scissors_exact_options(k, 0, 'c')


x = np.linspace(-4, 4, 100)
y = np.linspace(-4, 4, 100)

w1 = qt.wigner(sys.state, x, y)
cmap1 = qt.wigner_cmap(w1)

w2 = qt.wigner(statei, x, y)
cmap2 = qt.wigner_cmap(w2)

w3 = qt.wigner(sys2.state, x, y)
cmap3 = qt.wigner_cmap(w3)

w4 = qt.wigner(sys3.state, x, y)
cmap4 = qt.wigner_cmap(w4)

w5 = qt.wigner(sys4.state, x, y)
cmap5 = qt.wigner_cmap(w5)

fig, axes = plt.subplots(1, 5, figsize=(12,3))
im1 = axes[1].imshow(w1)
axes[1].set_title("Out")

norm = im1.norm

im0 = axes[0].imshow(w2, norm=norm)
axes[0].set_title("In")

im2 = axes[2].imshow(w3, norm=norm)
axes[2].set_title("a")

im3 = axes[3].imshow(w4, norm=norm)
axes[3].set_title("b")

im4 = axes[4].imshow(w5, norm=norm)
axes[4].set_title("c")

#fig.colorbar(im1)
fig.show()




fig, axes = plt.subplots(1, 5, figsize=(12,3))
qt.plot_fock_distribution(statei, fig=fig, ax=axes[0], title="In", unit_y_range=False);
qt.plot_fock_distribution(sys.state, fig=fig, ax=axes[1], title="Out", unit_y_range=False);
qt.plot_fock_distribution(sys2.state, fig=fig, ax=axes[2], title="a", unit_y_range=False);
qt.plot_fock_distribution(sys3.state, fig=fig, ax=axes[3], title="b", unit_y_range=False);
qt.plot_fock_distribution(sys4.state, fig=fig, ax=axes[4], title="c", unit_y_range=False);
fig.tight_layout()
plt.show()

