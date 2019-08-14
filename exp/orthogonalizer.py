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
N = 30
r = .6


## Initialize state
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r)

statei = sys.state.copy()
#a = qt.create(N) * qt.destroy(N)
a = qt.tensor(qt.create(N), qt.identity(N))
sys.ortho_oper(a)

sys.apply_orthogonalizer(b=0)

wplt.plot(sys.state, [-6, 6], 300)
print(sys.state.norm())

cm = sys.get_full_CM()

print(cm)
print(np.linalg.det(cm))
