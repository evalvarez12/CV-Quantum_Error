# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np
import measurements


N = 5
# Ref value for mu = 1.2 - EL = 2.7332
mu = 1.2

# Set initial state
sys = cv.System(N, 2)
r = np.arcsinh(np.sqrt(mu))
#sys.apply_SMS(r)
sys.apply_TMS(r)

# Make state no Gaussian
k = 0.4
sys.apply_photon_subtraction(k, 1)

psi0 = sys.state
if psi0.isket:
    rho0 = psi0 * psi0.dag()
    
# State with partial transpose
rho = qt.partial_transpose(rho0, [1, 0])

# Method 1 - Fastest
norm = rho.norm()
log_neg = np.log2(norm)

# Method 2    
rho2 = rho * rho.dag()
norm2 = (rho2.sqrtm()).tr()
log_neg2 = np.log2(norm2)

# Method 3 - for Gaussian states - best approx for small N?
sys.set_state(rho)
sys.set_quadratures_basis()
cm = sys.get_full_CM()
print(cm)
sym_eig = measurements.symplectic_eigenvalues(cm)
eig_ = np.min(sym_eig)
log_neg3 = max([0, -np.log2(eig_)])
#log_neg3 = -np.log2(eig_)


print(log_neg, log_neg2, log_neg3)