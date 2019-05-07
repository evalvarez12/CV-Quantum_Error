# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import tools
import cv_system as cv
import numpy as np
import qutip as qt
from scipy.linalg import block_diag


N = 2
sys = cv.System(N, Nmodes=2)
r = .6
sys.apply_SMD(r, 1)

ref_state = sys.state
print(ref_state)
   
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
r_aux = np.arcsinh(np.sqrt(0.5))
p = sys.apply_scissor(k, r_aux, 1)
print("P:", p)
print(sys.state)


