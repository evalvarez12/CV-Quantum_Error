# -*- coding: utf-8 -*-
"""
Testing PS in general over different inputs
Created on Thu Aug  8 10:14:13 2019

@author: Eduardo Villasenor 
"""

import src.cv_system as cv
import numpy as np
import qutip as qt
import wigner_plots as wplt

N = 20
r = 1
sys = cv.System(N, Nmodes=2)
sys.apply_TMS(r, [0, 1])
#sys.apply_SMS(.2)

#psi = sys.state
#a = qt.destroy(N)


cm = sys.get_full_CM()
print(cm)
#sys.apply_photon_subtraction(0.9)

#wplt.plot(a*psi, [-5,5])
#wplt.plot(sys.state, [-15,15], 300)
 

cm = sys.get_full_CM()
print(cm)

res = sys.homodyne_measurement()
print(res)

wplt.plot(sys.state, [-5,5], 300)
