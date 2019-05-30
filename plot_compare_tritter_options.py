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
import measurements


N = 10
sys = cv.System(N, 2)

eta = 0.01
m_aux = .01

k = 0.002
mu = 0.01

#sys.apply_SMD(r, 0)
r = np.arcsinh(np.sqrt(mu))
sys.apply_TMS(r)

sys.apply_loss_channel(eta, 1)

statei = sys.state
#statei = qt.basis(N, 1)   

r_aux = np.arcsinh(np.sqrt(m_aux))
#print(r_aux)
#p = sys.apply_scissors(k, r_aux, 1)
p = sys.apply_scissors_exact(k, 1)
#p11 = sys.apply_scissors(k, r_aux, 1)
#p = p*p11
#sys.apply_scissors_exact(k)




sys2 = cv.System(N, 1)

sys2.set_state(statei)
#pa = sys2.apply_scissors_options(k, r_aux, 1, 'a')
pa = sys2.apply_scissors_exact_options(k, 1, 'a')


sys3 = cv.System(N, 1)

sys3.set_state(statei)
#pb = sys3.apply_scissors_options(k, r_aux, 1, 'b')
pb = sys3.apply_scissors_exact_options(k, 1, 'b')


sys4 = cv.System(N, 1)

sys4.set_state(statei)
#pc = sys4.apply_scissors_options(k, r_aux, 1, 'c')
pc = sys4.apply_scissors_exact_options(k, 1, 'c')


sys_initial = cv.System(N, 2)
sys_initial.set_state(statei)


rci = measurements.CI(sys, [0])
rci_a = measurements.CI(sys2, [0])
rci_b = measurements.CI(sys3, [0])
rci_c = measurements.CI(sys4, [0])

rci_i = measurements.CI(sys_initial, [0])


#kr = measurements.key_rate(sys, 1, p)
#kra = measurements.key_rate(sys2, 1, pa)
#krb = measurements.key_rate(sys3, 1, pb)
#krc = measurements.key_rate(sys4, 1, pc)

print("rci:", rci)
print("rci_a:", rci_a)
print("rci_b:", rci_b)
print("rci_c:", rci_c)
print("rci_initial:", rci_i)



print("p:", p)
print("pa:", pa)
print("pb:", pb)
print("pc:", pc)

#print("KR:", kr)
#print("KRa:", kra)
#print("KRb:", krb)
#print("KRc:", krc)


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




#fig, axes = plt.subplots(1, 5, figsize=(12,3))
#qt.plot_fock_distribution(statei, fig=fig, ax=axes[0], title="In", unit_y_range=False);
#qt.plot_fock_distribution(sys.state, fig=fig, ax=axes[1], title="Out", unit_y_range=False);
#qt.plot_fock_distribution(sys2.state, fig=fig, ax=axes[2], title="a", unit_y_range=False);
#qt.plot_fock_distribution(sys3.state, fig=fig, ax=axes[3], title="b", unit_y_range=False);
#qt.plot_fock_distribution(sys4.state, fig=fig, ax=axes[4], title="c", unit_y_range=False);
#fig.tight_layout()
#plt.show()
#
