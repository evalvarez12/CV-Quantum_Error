# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np

# Define some initail parameters
N = 3

sys = cv.System(N, 1)

state = qt.tensor(qt.basis(N, 1), qt.basis(N, 1))
sys.add_state(state)
print(sys.state)
# Appl BS
theta = np.pi/4
posBS = [0, 1]
sys.apply_BS(theta, posBS)

print(sys.state)



sys2 = cv.System(N, 1)

sys2.add_state(state)
print(sys2.state)
posBS = [1, 0]
sys2.apply_BS(theta, posBS)
print(sys2.state)