# -*- coding: utf-8 -*-
"""
Fucntions recreate Braunstein CV-error correction

Created on Thu Feb 14 14:05:12 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import operatios as ops

# Define system paramters
N = 100
vacuum = qt.basis(N)

#Initial state
state = qt.displace(N, 1)*vacuum

# Define ancillas 
ancilla1 = qt.squeeze(0)*vacuum
ancilla2 = qt.squeeze(0)*vacuum


