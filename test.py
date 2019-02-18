# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 15:25:25 2019

@author: z5239621
"""
 
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import operations as ops

# Define some initail parameters
N = 100
eps = 1 + 2j
eps = eps/np.linalg.norm(eps)

S = qt.squeeze(N, eps)
D = qt.displace(N, eps)
vacuum = qt.basis(N)

state = S*vacuum
#state = D*vacuum


thetas = np.linspace(0, np.pi/2)

measurements = []
for i in thetas:
    measurements += [ops.S_homodyne(state, i)]


plt.plot(thetas, np.real(measurements))
plt.show()