# -*- coding: utf-8 -*-
"""
Test routines for quantum scissors

Created on Fri Mar  8 10:58:55 2019

@author: Eduardo Villasenor
"""
import sys
sys.path.append("..") 

import src.cv_system as cv
import numpy as np

N =10
k = .05
g = np.sqrt((1 - k)/k)
print('g:', g)



sys = cv.System(N, Nmodes=1, cm=False)
sys.apply_SMD(1)
print(sys.state)

sys.apply_scissors(k)
print(sys.state)


sys2 = cv.System(N, Nmodes=1, cm=False)
sys2.apply_SMD(1)


print(sys.state.dag() * sys2.state)