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


N = 5
sys = cv.System(N, )

r = .6
sys.apply_SMD(r, 0)

print(sys.state)
   
state = (qt.basis(N) + qt.basis(N, 1))/np.sqrt(2)
sys.add_state(state)



k = .01
sys.apply_scissor_exact(k, 1)
print(sys.state)



### SCISSOR NOT EXACT

N = 3
sys = cv.System(N, Nmodes=2)
r = .6
sys.apply_SMD(r, 1)

ref_state = sys.state
print(ref_state)
   
k = .01
m_aux = .1
r_aux = np.arcsinh(np.sqrt(m_aux))
p = sys.apply_scissor(k, r_aux, 1)
print("P:", p)
print(sys.state)

