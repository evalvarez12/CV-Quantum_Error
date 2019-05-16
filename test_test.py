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

statearr = sys.state.data.toarray()
inds = np.where(statearr != 0)[0]
print("indeces:", inds)
# Appl BS
theta = np.pi/4
posBS = [1, 2]
sys.apply_BS(theta, posBS)



print(sys.state)
statearr = sys.state.data.toarray()
inds = np.where(statearr != 0)[0]
print("indeces:", inds)




sys2 = cv.System(N, 1)

sys2.add_state(state)
print(sys2.state)
statearr2 = sys2.state.data.toarray()
inds2 = np.where(statearr2 != 0)[0]
print("indeces:", inds2)
posBS = [2, 1]
sys2.apply_BS(theta, posBS)
print(sys2.state)

statearr2 = sys2.state.data.toarray()
inds2 = np.where(statearr2 != 0)[0]
print("indeces:", inds2)