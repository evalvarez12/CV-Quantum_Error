# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np

# Define some initail parameters


N = 10
sys = cv.System(N, Nmodes=1)
r = .6
sys.apply_SMD(r)

ref_state = sys.state
k = .1
m_aux = .3
r_aux = np.arcsinh(np.sqrt(m_aux))
sys.apply_scissor(k, r_aux)

print(ref_state)
print(sys.state)
