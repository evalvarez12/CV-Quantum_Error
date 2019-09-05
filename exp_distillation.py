# -*- coding: utf-8 -*-
"""
File to test distillation of CAT_quasi_Bell states

Created on Wed Aug 21 14:13:06 2019

@author: Eduardo Villasenor
"""

import src.cv_system as cv
import numpy as np
import qutip as qt
import src.wigner_plots as wplt
N = 40
sys = cv.System(N, 0)
sys.add_CAT(3)

#sys.add_CAT_Bell(3, state='00')

qt.fock_distribution(qt.coherent(40, 3))



#wplt.plot(sys.state, [-7, 7], 200)