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



def start1(a, b, *argv):
    alpha, beta = argv
    Da = (alpha*a.dag() - alpha.conjugate()*a).expm()
    Db = (beta*b.dag() - beta.conjugate()*b).expm()
    return Da * Db


def start2(a, b):
    return a.dag()*b.dag()

def bs(a, b, theta):
    a_out = np.cos(theta)*a + 1j*np.sin(theta)*b
    b_out = np.cos(theta)*b + 1j*np.sin(theta)*a 
    return a_out, b_out


# Define some initial parameters
N = 5
t = .5
theta = np.arccos(np.sqrt(t))

alpha = .5
beta = .5

a = qt.tensor(qt.destroy(N), qt.identity(N))
b = qt.tensor(qt.identity(N), qt.destroy(N))
vacuum = qt.tensor(qt.basis(N), qt.basis(N))

Da = (alpha*a.dag() - alpha.conjugate()*a).expm()
Db = (beta*b.dag() - beta.conjugate()*b).expm()
D = Da * Db

U = operations.beam_splitter(N, theta)

a_out = np.cos(theta)*a + 1j*np.sin(theta)*b
b_out = np.cos(theta)*b + 1j*np.sin(theta)*a

Da_ref = (alpha*a_out.dag() - alpha.conjugate()*a_out).expm()
Db_ref = (beta*b_out.dag() - beta.conjugate()*b_out).expm()
D_ref = Da_ref * Db_ref

out = U * D * vacuum

out_ref = D_ref * vacuum

ins = D * vacuum 

#
#out = out.ptrace(1)
#out_ref = out_ref(1)

fig, axes = plt.subplots(1, 3, figsize=(12,3))
qt.plot_fock_distribution(D * vacuum, fig=fig, ax=axes[0], title="in", unit_y_range=False);
qt.plot_fock_distribution(out, fig=fig, ax=axes[1], title="og", unit_y_range=False);
qt.plot_fock_distribution(ins, fig=fig, ax=axes[2], title="ref", unit_y_range=False);
fig.tight_layout()
plt.show()


#N = 10
#sys = cv.System(N, 1)
#mu = 1.1
#r = np.arcsinh(np.sqrt(mu))
#sys.apply_SMD(r)
#
#ref_state = sys.state
#k = .1
#m_aux = .01
#r_aux = np.arcsinh(np.sqrt(m_aux))
##p = sys.apply_scissors(k, r_aux)
#p = sys.apply_scissors_options(k, r_aux, 0, 'c')
#
##sys.apply_photon_subtraction(k, 0)
#g = np.sqrt((1-k)/k)
#p_ref = np.exp((1-g**2)*np.linalg.norm(r)**2)/(1 + g**2)
#
#
#print(ref_state)
#print(sys.state)
#
#print("g:", g)
#print("p:", p)
#print("p_ref:", p_ref)
