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
   
k = .1
sys.apply_scissor_exact(k, 1)
print(sys.state)



state0 = qt.basis(N)
state1 = qt.basis(N, 1)
state = qt.tensor(state0, state1)
print(state)