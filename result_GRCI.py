# -*- coding: utf-8 -*-
"""
Replicating Gaussian RCI plots

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import numpy as np
import measurements

# Initial parameters
N = 10

# Sweep parameters
kappa = 0.005
mu = 0.05

# Fixed parameters
mu_aux = 0.01
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01


sys = cv.System(N, Nmodes=2, cm=False)
sys.save_state()
results = []
results2 = []
for kappa in np.linspace(0.0001, 0.01, 30):
    r2 = []
    r22 = []
    for mu in np.linspace(0.0001, .1, 30):
        
        r = np.arcsinh(np.sqrt(mu))
        sys.apply_TMS(r, pos=[0, 1])
        
        
        sys.apply_loss_channel(eta, 1)
        
        sys.apply_scissor(kappa, r_aux, 1)
        
        rci2 = measurements.RCI_simple(sys.state, [0])
        rci = measurements.RCI(sys, [0])
        print(mu, kappa, rci, rci2)
        r2 +=[rci]
        r22 += [rci2]
        
        sys.load_state()
    results += [r2]
    results2 += [r22]


filename = "data/rci_plot_1NLA"

results = np.array(results)
np.save(filename, results)


filename2 = "data/rci_plot_1NLA2"

results2 = np.array(results2)
np.save(filename2, results2)




#
##results = []
##for kappa in np.linspace(0.001, 0.01, 20):
##    r2 = []
##    for mu in np.linspace(0.001, .1, 20):
#
#
#vacuum = qt.basis(N) * qt.basis(N).dag()
#vacuum = qt.tensor(vacuum, vacuum)
#
#s = np.arcsinh(np.sqrt(mu))
#S =  ops.tmsqueeze(N, s)
#
#rho_in = S * vacuum * S.dag()
#
#rho = bs.loss_channel_applyU(rho_in, 0, eta)
#
#rho = scissor.NLA(rho, kappa, mu_aux)
#
#
#rci2 = ops.RCI(rho, 1)
#
#
#q = (qt.destroy(N) + qt.destroy(N).dag())/np.sqrt(2)
#p = (qt.destroy(N) - qt.destroy(N).dag())/(1j*np.sqrt(2))
#
#
##  symplectic eigenvalues over A
#cm_basis = [q, p]
#
#rhoA = rho.ptrace(0)
#cm = 2 * qt.covariance_matrix(cm_basis, rhoA, symmetrized=True).astype(float)
#
## Chack if symetric
#print(cm == cm.transpose())
#
## Chack if positive-definite
#print(is_pos_def(cm))
#
#omega = np.array([[0, 1], [-1, 0]])
#
#symplectic_eigA =np.linalg.eigvals(1j * np.dot(omega, cm))
#
#
#
## symplectic eigenvalues over AB
#q1 = tools.tensor(q, N, 0, 2)
#q2 = tools.tensor(q, N, 1, 2)
#p1 = tools.tensor(p, N, 0, 2)
#p2 = tools.tensor(p, N, 1, 2)
#
#cm_basis = [q1, p1, q2, p2]
#
#cm = 2 * qt.covariance_matrix(cm_basis, rho, symmetrized=True).astype(float)
#
## Chack if symetric
#print(cm == cm.transpose())
#
## Chack if positive-definite
#print(is_pos_def(cm))
#
#omega = np.array([[0, 1], [-1, 0]])
#omega = la.block_diag(omega, omega)
#
#symplectic_eigAB =np.linalg.eigvals(1j * np.dot(omega, cm))


# Calculate the RCI
#rci = g(max(symplectic_eigA).real) -  g(max(symplectic_eigAB).real)


#        r2 +=[rci]
#    results += [r2]
#
#
#filename = "rci_plot_1NLA"
#
#results = np.array(results)
#np.save(filename, results)






###### Sweep over paramter space
#results = []
#
#for kappa in np.linspace(0.001, 0.01, 20):
#    r2 = []
#    for mu in np.linspace(0.001, .1, 20):
#        
#        s = np.arcsinh(np.sqrt(mu))
#        S =  ops.tmsqueeze(N, s)
#        
#        rho_in = S * vacuum * S.dag()
#        
#        rho = bs.loss_channel_applyU(rho_in, 0, eta)
#        
#        rho = scissor.NLA(rho, kappa, mu_aux)
#        
#        
#        rci = ops.RCI(rho, 1)
#
#        r2 +=[rci]
#    results += [r2]
#
#
#
#filename = "rci_plot_1NLA"
#
#results = np.array(results)
#np.save(filename, results)
