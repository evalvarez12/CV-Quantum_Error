# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 10:32:06 2020

@author: z5239621
"""

import numpy as np
import qutip as qt
import src.cv_system as cv
import src.wigner_plots as wp

N = 20
r = 1

sys = cv.System(N, 0)
sys.add_fock(1)
#sys.apply_SMS(r)
#sys.apply_SMD(r)

#sys.apply_photon_subtraction(.8)
cm = sys.get_full_CM()
print(cm)


wp.plot(sys.state, [-5, 5])