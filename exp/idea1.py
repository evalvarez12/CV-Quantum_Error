# -*- coding: utf-8 -*-
"""
Test the idea of an entangled photon to a qubit interacting with a CV-state

Created on Wed May  8 15:07:58 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import src.cv_system as cv
import wigner_plots as wplt

sys.path.append("..") 



N = 20
sys = cv.System(N, 1)

r = .5
sys.apply_SMS(r)
sys.apply_photon_catalysis(1, .9)
#print(statei)
  
wplt.plot(sys.state, [-4,4])