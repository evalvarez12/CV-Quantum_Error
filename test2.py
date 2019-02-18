# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:36:08 2019

@author: z5239621
"""

import qutip as qt
import numpy as np
import operations as ops

# Define some initail parameters
N = 3

#S = qt.squeeze(N, eps)
#D = qt.displace(N, 1)
vacuum = qt.basis(N)

state_zero = qt.tensor(vacuum, vacuum)

a = qt.destroy(N)
b = qt.destroy(N)

a,b = ops.beam_splitter([a, b], np.pi/4)

state = a.dag()*b.dag()*state_zero

print(state)