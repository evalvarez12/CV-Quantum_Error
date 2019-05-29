# -*- coding: utf-8 -*-
"""
Files to the test photon subtraction operation

Created on Thu Apr  4 15:56:52 2019

@author: Eduardo Villasenor 
"""
import cv_system as cv
import numpy as np
import qutip as qt
from scipy.special import comb

N = 50
T = 0.9
mean_photon_number = 5

r = np.arcsinh(np.sqrt(mean_photon_number))
print("Squeezing paramter:", r)
# TMSV after photon subtraction
def alpha_n(r, n):
    return (np.tanh(r)**n)/np.cosh(r)

def r_nkT(n, k , t):
    return np.sqrt(comb(n, k) * t**(n-k) * (1 - t)**k )

psi_theory = qt.tensor(qt.basis(N), qt.basis(N)) * 0

for i in range(N):
#    psi_theory += alpha_n(i) * r_nkT(i, 1, T) * qt.tensor(qt.basis(N, i), qt.basis(N, i))
    psi_theory += alpha_n(r, i) * qt.tensor(qt.basis(N, i), qt.basis(N, i))
    print(i, alpha_n(r, i))


norm = psi_theory.norm()
psi_theory = psi_theory/norm
print('Norm:', norm)


sys = cv.System(N, Nmodes=2)
sys.apply_TMS(mean_photon_number, [0, 1])


psi_simulation = sys.state


a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))

n_op = a.dag()*a 

print("Mean photon number theory:", qt.expect(n_op, psi_theory))
print("Mean photon number simulation:", qt.expect(n_op, psi_simulation))



fig, axes = plt.subplots(1, 2, figsize=(12,3))
qt.plot_fock_distribution(psi_theory, fig=fig, ax=axes[0], title="Theory", unit_y_range=False);
qt.plot_fock_distribution(psi_simulation, fig=fig, ax=axes[1], title="Simulated", unit_y_range=False);
#fig.tight_layout()

#fig, axes = plt.subplots()
#qt.plot_fock_distribution(psi_theory, fig=fig, ax=axes, title="Theory", unit_y_range=True);

plt.show()