# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import operations as ops
import beam_splitter as bs

N = 10

vacuum = qt.basis(N, 0)

meanpn = 5
r = np.arcsinh(np.sqrt(meanpn))

stateAB = ops.tmsqueeze(N, r) * qt.tensor(vacuum, vacuum)

stateABC = qt.tensor(stateAB, vacuum)

# Transissivty of photon subtraction BS
tps = 0.9
stateABC = bs.beam_splitter_Uoperator(N, tps, [1, 2], 3) * stateABC

stateAB = ops.collapse_photon_number(N, stateABC, n=1, pos=2, Nmodes=3)




