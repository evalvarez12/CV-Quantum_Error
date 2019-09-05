# -*- coding: utf-8 -*-
"""
Testing shared information calculation method
Created on Tue Aug 27 17:05:55 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import src.cv_system as cv
import src.measurements as measurements


N = 20
sys = cv.System(N, Nmodes=2, cm=False)
r = .5
sys.apply_TMS(r, [0, 1])

sys.set_quadratures_basis()
Va = sys.get_CM_entry([0, 0])
Vb = sys.get_CM_entry([2, 2])

Cab = sys.get_CM_entry([0, 2])

I = measurements.I(Va, Vb, Cab)

I2 = measurements.I_altern(sys.state, [1, 0])

print(I, I2)