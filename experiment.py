# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np
import operations
import matplotlib.pyplot as plt
import hamiltonians as ham
import symplectic as sym
import theory
import measurements


N = 10
mu = 2

sys = cv.System(N, 2)
r = np.arcsinh(np.sqrt(mu))
#sys.apply_SMS(r)
sys.apply_TMS(r)

psi0 = sys.state

if psi0.isket:
    rho0 = psi0 * psi0.dag()
    
rho = qt.partial_transpose(rho0, [0, 1])

norm = rho.norm()
log_neg = np.log2(norm)
    
rho2 = rho * rho.dag()
norm2 = (rho2.sqrtm()).tr()
log_neg2 = np.log2(norm2)


print(norm, norm2)

sys.set_state(rho)
sys.set_quadratures_basis()
cm = sys.get_full_CM()
print(cm)
sym_eig = measurements.symplectic_eigenvalues(cm)
eig_ = np.min(sym_eig)
log_neg3 = max([0, -np.log2(eig_)])
#log_neg3 = -np.log2(eig_)


print(log_neg, log_neg2, log_neg3)