# -*- coding: utf-8 -*-
"""
Create results for the GLOBECOMM paper
Created on Mon Mar  8 11:32:18 2021

@author: z5239621
"""


import TMSV_sq as tmsv
import qsc as qs_t
import qsc_rec as qs_r
import sb_sq as sb
import scipy.io
import numpy as np
import matplotlib.pyplot as plt


#T = np.random.rand()
# T = 1
#s = np.random.rand()
eta = np.sqrt(10**(-1/10))
g = 1 / eta
eps = 0
V = np.random.rand()*5 + 1
d = np.pi * np.random.rand()*2

T = 0.9
s = 0.7
print(T, s)

# V=1.909335096180439
#print('TMSV', tmsv.fidelity(V, T, eps, g, eta, s)) 
#print('SB', sb.fidelity(V, T, d, eps, eta, g, s))

print('-----------------------------------------------------------')

print('TMSV', tmsv.opt_fidelity(T, eps, eta, s))
print('SB', sb.opt_fidelity(T, eps, eta, s))


print('SB', sb.opt_fidelity_avg(T, eps, eta, s, .1))
print('TMSV', tmsv.opt_fidelity_avg(T, eps, eta, s, .1))

#print('TMSV ---', tmsv.fidelity(18.91863554, T, eps, g, eta, 0.7341))

############################# FIXED LOSS FIDELITIES PLOT
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# T = np.linspace(0.01, 0.99, 15)
# # T = [0.74, 0.81, 0.9, 1]
# #T = np.logspace(-3, 0, 15)
# eps = 0.05
# eta = np.sqrt(10**(-1/10))
# s = [.575]
# # eta = 1
# #eta = 0.7
#
# f_tmsv = []
# f_qs_t = []
# f_sb = []
#
# for sig_ind in range(len(s)):
#     f_tmsv_i = []
#     f_qs_t_i = []
#     f_sb_i = []
#
#     for it in T:
#         f_tmsv_i += [tmsv.opt_fidelity(it, eps, eta, s[sig_ind])]
#         f_qs_t_i += [qs_t.opt_avg_fidelity(it, eps, eta, s[sig_ind])]
# #        f_qs_r_i += [qs_r.opt_avg_fidelity(it, eps, eta, sigmaT[sig_ind])]
#         f_sb_i += [sb.opt_fidelity(it, eps, eta, s[sig_ind])]
#
#
#     f_tmsv += [f_tmsv_i]
#     f_qs_t += [f_qs_t_i]
#     f_sb += [f_sb_i]
#
#
#
#
# ax.plot(T, f_tmsv[0], label=r'TMSV $\sigma= $' + str(np.round(np.abs(s[0]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
# # ax.plot(T, f_tmsv[1], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[1]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
# # ax.plot(T, f_tmsv[2], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[2]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
# # ax.plot(T, f_tmsv[3], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
#
# #ax.plot(T, f_qs_t[0], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[0]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
# #ax.plot(T, f_qs_t[1], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[1]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
# #ax.plot(T, f_qs_t[2], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[2]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
# #ax.plot(T, f_qs_t[3], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#
# ax.plot(T, f_sb[0], label=r'SB_r $\sigma= $' + str(np.round(np.abs(s[0]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
# # ax.plot(T, f_sb[1], label=r'SB_r $\sigma= $' + str(np.round(np.abs(s[1]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
# # ax.plot(T, f_sb[2], label=r'SB_r $\sigma= $' + str(np.round(np.abs(s[2]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
# # ax.plot(T, f_sb[3], label=r'SB_r $\sigma= $' + str(np.round(np.abs(s[3]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
#
#
# clasical = np.ones_like(f_tmsv[0]) * .5
# plt.plot(T, clasical, 'r--')
#
# ax.set_xlabel(r'T')
# ax.set_ylabel(r'$\bar{\mathcal{F}} $')
# # plt.ylim([-0.005, 0.2])
#
# #ax.set_xlim(-0.1, 13)
# #ax.set_ylim(.49, 1)
# #ax.set_xscale('log')
# ax.grid()
# ax.legend()
# #ax.set_xscale('log')
# ### To MATLAB
# data = [T, f_tmsv, f_qs_t, f_sb]
# #scipy.io.savemat('figs_noG/fixed_py_eta0.7.mat', {'data':data})
#
# plt.show()
