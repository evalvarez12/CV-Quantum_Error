# -*- coding: utf-8 -*-
"""
Code to test CV system
Created on Mon Apr 15 10:47:50 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import tools
import matplotlib.pyplot as plt
# Parameters
N = 40
mpn = 1
t = .5

r = np.arcsinh(np.sqrt(mpn))
print("Squeezing:", r)

## Initialize state
sys = cv.System(N, Nmodes=2)
sys.apply_TMS(mpn, [0, 1])
# BAD TMSV
#sys.replace_current_state_w_bad_TMSV(mean_photon_number)


sys.set_quadratures_basis()
CM = sys.get_full_CM()
print("CM:", CM)

print("MG simple")
Va = sys.get_simple_CM_V(0)
Vb = sys.get_simple_CM_V(1)
Cab = sys.get_simple_CM_C([0, 1])
print("Va:", Va)
print("Vb:", Vb)
print("Cab:", Cab)

print("Theory comparisson")
print("V:", 2*mpn+1)
print("C:", 2*np.sqrt(mpn*(mpn+1)))
#sys.apply_BS(t, [0, 1])

#Vc = sys.get_simple_CM_V(2)

#Cbc = sys.get_simple_CM_C([1, 2])
#Cac = sys.get_simple_CM_C([0, 2])


statea = sys.state


q1 = sys.quad_basis[0]
p1 = sys.quad_basis[1]
q2 = sys.quad_basis[2]
p2 = sys.quad_basis[3]

print("Computed:", (qt.expect(q2*q1 + q1*q2, sys.state) - 2*qt.expect(q1, sys.state)*qt.expect(q2, sys.state)))

print("-----------------> Apply BS")

sys.apply_BS(t, [0, 1])

sys.set_quadratures_basis()
CM = sys.get_full_CM()
print("CM:", CM)

print("MG simple")
Va = sys.get_simple_CM_V(0)
Vb = sys.get_simple_CM_V(1)
Cab = sys.get_simple_CM_C([0, 1])
print("Va:", Va)
print("Vb:", Vb)
print("Cab:", Cab)

print("Theory comparisson")
print("V1:", np.exp(2*r))
print("V2:", np.exp(-2*r))


print("-----------------> Sys 2")
sys2 = cv.System(N, Nmodes=2)
#z = mpn*np.exo(1j*np.pi/4)
sys2.apply_SMS(mpn, 0)
sys2.apply_SMS(mpn, 1)

sys2.set_quadratures_basis()
CM2 = sys2.get_full_CM()
print("CM:", CM2)


q1 = sys2.quad_basis[0]
p1 = sys2.quad_basis[1]
q2 = sys2.quad_basis[2]
p2 = sys2.quad_basis[3]

print("Computed:", (qt.expect(p2*p2 + p2*p2, 1j*sys2.state) - 2*qt.expect(p2, 1j*sys2.state)*qt.expect(p2, 1j*sys2.state)))


sys2.apply_BS(t, [0, 1])
sys2.set_quadratures_basis()
CM2 = sys2.get_full_CM()
print("CM:", CM2)

#
#fig, axes = plt.subplots(1, 3, figsize=(12,3))
#qt.plot_fock_distribution(statea, fig=fig, ax=axes[0], title="A", unit_y_range=False);
#qt.plot_fock_distribution(sys.state, fig=fig, ax=axes[1], title="B", unit_y_range=False);
#qt.plot_fock_distribution(sys2.state, fig=fig, ax=axes[2], title="C", unit_y_range=False);
#
##fig.tight_layout()
#plt.show()
#









#########################################################################################
#sys.set_quadratures_basis()
#CM = sys.get_full_CM()
#print("CM:", CM)



#Va = sys.get_simple_CM_V(0)
#Vb = sys.get_simple_CM_V(1)
##Vc = sys.get_simple_CM_V(2)
#
#Cab = sys.get_simple_CM_C([0, 1])
##Cbc = sys.get_simple_CM_C([1, 2])
##Cac = sys.get_simple_CM_C([0, 2])
#
#
#print("Va:", Va)
#print("Vb:", Vb)
##print("Vc:", Vc)
#
#print("Cab:", Cab)


#sys = cv.System(N, Nmodes=2)
#sys.apply_SMS(mpn, 0)
#
#sys.set_quadratures_basis()
#CM = sys.get_full_CM()
#print("CM:", CM)
#a = qt.destroy(N)
#a1 = tools.tensor(N, a, 1, 3)
#a2 = tools.tensor(N, a, 2, 3)
#c = sys.state.dag() * (a1*a2 + a1.dag()*a2.dag()) * sys.state
#print(c.norm())
#
#sys.collapse_fock_state(0, 2)
#
#Va = sys.get_simple_CM_V(0)
#Vb = sys.get_simple_CM_V(1)
#
#Cab = sys.get_simple_CM_C([0, 1])
#
#
#print("Va:", Va)
#print("Vb:", Vb)
#
#
#print("Cab:", Cab)
#
#sys.set_quadratures_basis()
#CM = sys.get_full_CM()
#print("CM:", CM)

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