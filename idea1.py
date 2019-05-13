# -*- coding: utf-8 -*-
"""
Test the idea of an entangled photon to a qubit interacting with a CV-state

Created on Wed May  8 15:07:58 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import tools
import matplotlib.pyplot as plt


#N = 5
#sys = cv.System(N, 1)
#
#r = .6
#sys.apply_SMD(r, 0)
#statei = sys.state
#print(statei)
#   
#
#k = 1
#s0 = qt.tensor(qt.basis(N), qt.basis(N))
#s1 = qt.tensor(qt.basis(N, 1), qt.basis(N, 1))
#
##state = (s0 + k*s1)
##state = state/state.norm()
#state = (qt.basis(N) + qt.basis(N, 1))/np.sqrt(2)
#sys.add_state(state)
#
#print(state)
#
#print(sys.state)
#
#t = 0.5
#theta = np.arccos(np.sqrt(t))
#sys.apply_BS(theta, pos=[0, 1])
#
#print(sys.state)
#
#sys.collapse_fock_state(0, 1)
#
#print(sys.state)
#print(statei)




N = 20
sys = cv.System(N, 1)

r = .5
sys.apply_SMD(r, 0)
statei = sys.state
print(statei)
   

k = 1
s0 = qt.tensor(qt.basis(N), qt.basis(N))
s1 = qt.tensor(qt.basis(N, 1), qt.basis(N, 1))

state = (s0 + k*s1)
state = state/state.norm()
sys.add_state(state)

print(state)

print(sys.state)

t = 0.5
theta = np.arccos(np.sqrt(t))
sys.apply_BS(theta, pos=[0, 1])

focks = [1, 1]

sys.collapse_fock_state(focks[0], 1)
print(sys.state)

sys.collapse_fock_state(focks[1], 1)

print(statei)
print(sys.state)

x = np.linspace(-3, 5, 100)
y = np.linspace(-4, 4, 100)

w1 = qt.wigner(sys.state, x, y)
cmap1 = qt.wigner_cmap(w1)

w2 = qt.wigner(statei, x, y)
cmap2 = qt.wigner_cmap(w2)

fig, axes = plt.subplots(1, 2, figsize=(12,3))
im1 = axes[1].imshow(w1)
axes[1].set_title("Out")

norm = im1.norm

im0 = axes[0].imshow(w2, norm=norm)
axes[0].set_title("In")
fig.colorbar(im1)
fig.show()



#plt.subplot(1, 4, 1)
#plt.imshow(w1, cmap1)
#plt.colorbar()
#plt.subplot(1, 4, 2)
#plt.imshow(w2, cmap2)
#plt.colorbar()

#plt.show()

