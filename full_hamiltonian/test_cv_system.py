# -*- coding: utf-8 -*-
"""
Code to test CV system
Created on Mon Apr 15 10:47:50 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import tools

# Parameters
N = 20
mpn = 1.3
t = .5


## Initialize state
sys = cv.System(N, Nmodes=2)
sys.apply_TMS(mpn, [0, 1])
# BAD TMSV
#sys.replace_current_state_w_bad_TMSV(mean_photon_number)


sys.add_vacuum()
sys.apply_BS(t, [1, 2])

Va = sys.get_simple_CM_V(0)
Vb = sys.get_simple_CM_V(1)
Vc = sys.get_simple_CM_V(2)

Cab = sys.get_simple_CM_C([0, 1])
Cbc = sys.get_simple_CM_C([1, 2])
Cac = sys.get_simple_CM_C([0, 2])


print("Va:", Va)
print("Vb:", Vb)
print("Vc:", Vc)

print("Cab:", Cab)
print("Cbc:", Cbc)
print("Cac:", Cac)

CM = sys.get_full_CM()
print("CM:", CM)


a = qt.destroy(N)
a1 = tools.tensor(N, a, 1, 3)
a2 = tools.tensor(N, a, 2, 3)
c = sys.state.dag() * (a1*a2 + a1.dag()*a2.dag()) * sys.state
print(c.norm())

#a1 = tools.tensor(N, qt.create(N), 0, 3)
#a2 = tools.tensor(N, qt.create(N), 1, 3)
#a3 = tools.tensor(N, qt.create(N), 2, 3)
#
#permute_list = [0, 2, 1]
#a3 = a3.permute(permute_list)
#
#sys.state = a2 * sys.state
#
#Va = sys.get_simple_CM_V(0)
#Vb = sys.get_simple_CM_V(1)
#Vc = sys.get_simple_CM_V(2)
#
#print("------------------------------")
#print("Va:", Va)
#print("Vb:", Vb)
#print("Vc:", Vc)