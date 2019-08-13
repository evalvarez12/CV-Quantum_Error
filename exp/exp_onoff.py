# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np
import operations
import matplotlib.pyplot as plt
import hamiltonians as ham
import symplectic as sym
import theory

N = 6

a = np.random.rand(N)
a = a/np.linalg.norm(a)


statei = qt.basis(N)*0
for i in range(N):
    statei += a[i]*qt.basis(N, i)

sys = cv.System(N)
sys.set_state(statei)
sys.add_TMSV(.1)


p = sys.collapse_ON_OFF(measurement=0, pos=2)
print(statei)

print("p:", p, np.sum(a[1:]**2))
print(sys.state)

if not sys.state.isket:
    print("purity:", (sys.state* sys.state).tr())