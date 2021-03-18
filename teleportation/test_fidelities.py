# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:15:00 2021

@author: z5239621
"""

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import TMSV as tmsv
import qsc as qs
import sb as sb


V = np.random.rand() + 1
T = np.random.rand() * 1
tsc = np.random.rand() 
#T = 1
#V = 1.00001
eps = 0.05
alpha = 5
sigma = 5
eta = 1
g = 1
print('-------------Random')
print('pars:', np.round([V, T, tsc], 3))
print('F_TMSV:', tmsv.fidelity(V, T, eps, eta, g, alpha))
print('F_TMSV_avg:', tmsv.fidelity_alphabet(V, T, eps, eta, g, sigma))
print('F_QS:', qs.fidelity(V, T, tsc, eps, eta, g, alpha))
print('F_QS_avg:', qs.avg_fidelity(V, T, tsc, eps, eta, g, sigma))
print('F_SB:', sb.squeezed_bell(V, T, tsc, eps, eta, g, alpha))
print('F_SB_avg:', sb.squeezed_bell_avg(V, T, tsc, eps, eta, g, sigma))

print('-------------------')
#print('F_opt', opt_fidelity(T, eps, 1, alpha))

eps = 0.0
alpha = 5

sigma = 1
eta = 1
g = 1
T = 0.76
tsc = 0.7
V = 1.8
alpha = 1
gsc = np.sqrt((1 - tsc)/tsc)
#print('gsc:', gsc)
#print('p_succ:', 2 / (V+1) - (gsc*2/(V+1))**2 * (1 / (1 + gsc**2)))
#
#print('F:', fidelity(V, T, tsc, eps, eta, g, alpha))
#print('F_opt', opt_fidelity(T, eps, eta, sigma))
print('F_TMSV:', tmsv.fidelity(V, T, eps, eta, g, alpha))
print('F_avg:', qs.avg_fidelity(V, T, tsc, eps, eta, g, sigma))
print('F_SB_avg:', sb.squeezed_bell_avg(V, T, tsc, eps, eta, g, sigma))

#print('F_avg_TMSV:', tmsv.opt_fidelity_alphabet_gopt(T, eps, eta, sigma))
#print('F_avg_opt', opt_avg_fidelity(T, eps, eta, sigma))
#print('F hand:', fidelity_pars([1.6, .999, .541], T, eps, 1, alpha))