# -*- coding: utf-8 -*-
"""
Functions to test orthogonalizer

Created on Tue Aug 13 11:32:13 2019

@author: Eduardo Villasenor
"""

import src.cv_system as cv
#import plots.wigner_plots as wplt
import wigner_plots as wplt
import qutip as qt


# Parameters
N = 20
r = .7


## Initialize state
sys = cv.System(N, Nmodes=1, cm=False)
sys.apply_SMS(r)

statei = sys.state.copy()


#a = qt.destroy(N) 
##a = qt.tensor(qt.create(N), qt.identity(N))
#sys.ortho_oper(a)
#sys.apply_orthogonalizer(b=.4)


#sys.apply_photon_catalysis(0.5, 1)
#sys.apply_photon_subtraction(0.9)


#x_res = sys.homodyne_measurement()
#print("X:", x_res)

wplt.plot(sys.state, [-6, 6], 200)
print(sys.state.norm())
print(statei.dag() * sys.state)

#cm = sys.get_full_CM()
#
#print(cm)
#print(np.linalg.det(cm))
