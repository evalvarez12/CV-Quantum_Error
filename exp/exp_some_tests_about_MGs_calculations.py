# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import src.cv_system as cv
import qutip as qt
import numpy as np



N = 30
mu = 2

sys = cv.System(N, 2)


r = 1
sys.apply_TMS(r)
#sys.apply_SMS(r, 1)
#sys.apply_SMD(r, 1)


#state = qt.basis(N, 1)
#sys.add_state(state)

for i in range(1):
    sys.apply_photon_subtraction(0.95, 0)

#sys.apply_loss_channel(.1, 1)


#sys.apply_scissors_exact(0.05, 0)
#sys.apply_photon_catalysis(4, 0.9, 0)
#sys.apply_photon_catalysis(3, 0.9, 0)
#sys.apply_photon_catalysis(2, 0.9, 0)
#sys.apply_photon_catalysis(1, 0.9, 0)


sys.set_quadratures_basis()

sys.set_full_CM()

print(sys.cm)
print(np.linalg.det(sys.cm))
c_aa = sys.get_simple_CM_V(0)
c_bb = sys.get_simple_CM_V(1)
c_ab = sys.get_simple_CM_C([0, 1])

print(c_aa)
print(c_bb)
print(c_ab)

#print("---------------------------------------")
#sys.apply_BS(np.pi/6, [0, 1])
#print(sys.cm)
#
#cm = sys.get_full_CM()
#print(cm)
#
#c_aa = sys.get_simple_CM_V(0)
#c_bb = sys.get_simple_CM_V(1)
#c_ab = sys.get_simple_CM_C([0, 1])
#
#print(c_aa)
#print(c_bb)
#print(c_ab)