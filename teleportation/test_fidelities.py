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


V = 1
#T = np.logspace(-3, 0, 15)
T = 0.9999
eps = 0.005
alpha = 5
eta = .5
g = 2 
d = 0.66934996
print(sb.fidelity(V, T, d, eps, eta, g, alpha))
print(sb.opt_fidelity_r(V, T, eps, eta, alpha))


if True:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    V = np.linspace(1, 10, 15)
    #T = np.logspace(-3, 0, 15)
    T = [0.1, 0.4, 0.7, .99]
    eps = 0.005
    eta = 10**(-1/10)
    alpha = 5
    
    f_tmsv = []
    f_sb = []
    
    g = 1
    for iT in T:
        f_tmsv_i = []
        f_sb_i = []
        
        for iV in V:
            f_tmsv_i += [tmsv.opt_fidelity_r(iV, iT, eps, eta, alpha)]
            f_sb_i += [sb.opt_fidelity_r(iV, iT, eps, eta, alpha)]
    #        f_tmsv_i += [tmsv.opt_avg_fidelity_r(iV, iT, eps, eta, alpha)]
    #        f_sb_i += [sb.opt_avg_fidelity_r(iV, iT, eps, eta, alpha)]
    
        
        f_tmsv += [f_tmsv_i]
        f_sb += [f_sb_i]
    
        
        
        
    ax.plot(V, f_tmsv[0], label=r'TMSV T= ' + str(np.round(np.abs(T[0]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
    ax.plot(V, f_tmsv[1], label=r'TMSV T= ' + str(np.round(np.abs(T[1]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
    ax.plot(V, f_tmsv[2], label=r'TMSV T= ' + str(np.round(np.abs(T[2]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
    ax.plot(V, f_tmsv[3], label=r'TMSV T= ' + str(np.round(np.abs(T[3]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
    
    ax.plot(V, f_sb[0], label=r'SB T= ' + str(np.round(np.abs(T[0]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
    ax.plot(V, f_sb[1], label=r'SB T= ' + str(np.round(np.abs(T[1]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
    ax.plot(V, f_sb[2], label=r'SB T= ' + str(np.round(np.abs(T[2]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
    ax.plot(V, f_sb[3], label=r'SB T= ' + str(np.round(np.abs(T[3]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
    
    
    clasical = np.ones_like(f_tmsv[0]) * .5
    plt.plot(V, clasical, 'r--')
    
    ax.set_xlabel(r'V')
    ax.set_ylabel(r'$\bar{\mathcal{F}} $')
    # plt.ylim([-0.005, 0.2])
    
    #ax.set_xlim(-0.1, 13)
    #ax.set_ylim(.49, 1)
    #ax.set_xscale('log')
    ax.grid()
    ax.legend()
    
    plt.show()
    












####################################################################
####################################################################
#V = np.random.rand() + 1
#T = np.random.rand() * 1
#tsc = np.random.rand() 
##T = 1
##V = 1.00001
#eps = 0.05
#alpha = 5
#sigma = 5
#eta = 1
#g = 1
#print('-------------Random')
#print('pars:', np.round([V, T, tsc], 3))
#print('F_TMSV:', tmsv.fidelity(V, T, eps, eta, g, alpha))
#print('F_TMSV_avg:', tmsv.fidelity_alphabet(V, T, eps, eta, g, sigma))
#print('F_QS:', qs.fidelity(V, T, tsc, eps, eta, g, alpha))
#print('F_QS_avg:', qs.avg_fidelity(V, T, tsc, eps, eta, g, sigma))
#print('F_SB:', sb.squeezed_bell(V, T, tsc, eps, eta, g, alpha))
#print('F_SB_avg:', sb.squeezed_bell_avg(V, T, tsc, eps, eta, g, sigma))
#
#print('-------------------')
##print('F_opt', opt_fidelity(T, eps, 1, alpha))
#
#eps = 0.0
#alpha = 5
#
#sigma = 1
#eta = 1
#g = 1
#T = 0.76
#tsc = 0.7
#V = 1.8
#alpha = 1
#gsc = np.sqrt((1 - tsc)/tsc)
##print('gsc:', gsc)
##print('p_succ:', 2 / (V+1) - (gsc*2/(V+1))**2 * (1 / (1 + gsc**2)))
##
##print('F:', fidelity(V, T, tsc, eps, eta, g, alpha))
##print('F_opt', opt_fidelity(T, eps, eta, sigma))
#print('F_TMSV:', tmsv.fidelity(V, T, eps, eta, g, alpha))
#print('F_avg:', qs.avg_fidelity(V, T, tsc, eps, eta, g, sigma))
#print('F_SB_avg:', sb.squeezed_bell_avg(V, T, tsc, eps, eta, g, sigma))
#
##print('F_avg_TMSV:', tmsv.opt_fidelity_alphabet_gopt(T, eps, eta, sigma))
##print('F_avg_opt', opt_avg_fidelity(T, eps, eta, sigma))
##print('F hand:', fidelity_pars([1.6, .999, .541], T, eps, 1, alpha))