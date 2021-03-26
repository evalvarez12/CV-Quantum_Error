# -*- coding: utf-8 -*-
"""
Create results for the GLOBECOMM paper
Created on Mon Mar  8 11:32:18 2021

@author: z5239621
"""


import TMSV as tmsv
import qsc as qs_t
import qsc_rec as qs_r
import sb as sb
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

# plt.close('all')

############################## FIXED LOSS FIDELITIES PLOT - NON AVERAGE
#fig = plt.figure()
#ax = fig.add_subplot(111)
#
#T = np.linspace(0.01, 1)
#eps = 0.02
#eta = 10**(-1/10)
#alpha = [1, 3, 8, 10]
#
#f_tmsv = []
#f_qs_t = []
#f_qs_r = []
#
#for alp_ind in range(len(alpha)):
#    f_tmsv_i = []
#    f_qs_t_i = []
#    f_qs_r_i = []
#
#    for it in T:
#        f_tmsv_i += [tmsv.opt_fidelity(it, eps, eta, alpha[alp_ind])]
#        f_qs_t_i += [qs_t.opt_fidelity(it, eps, eta, alpha[alp_ind])]
##        f_qs_r_i += [qs_r.opt_avg_fidelity(it, eps, eta, sigmaT[sig_ind])]
#
#    f_tmsv += [f_tmsv_i]
#    f_qs_t += [f_qs_t_i]
#    f_qs_r += [f_qs_r_i]
#
#
#
#
#ax.plot(T, f_tmsv[0], label=r'TMSV $\alpha= $' + str(np.round(np.abs(alpha[0]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
#ax.plot(T, f_tmsv[1], label=r'TMSV $\alpha= $' + str(np.round(np.abs(alpha[1]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
#ax.plot(T, f_tmsv[2], label=r'TMSV $\alpha= $' + str(np.round(np.abs(alpha[2]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
#ax.plot(T, f_tmsv[3], label=r'TMSV $\alpha= $' + str(np.round(np.abs(alpha[3]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
#
#ax.plot(T, f_qs_t[0], label=r'QS_t $\alpha= $' + str(np.round(np.abs(alpha[0]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[1], label=r'QS_t $\alpha= $' + str(np.round(np.abs(alpha[1]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[2], label=r'QS_t $\alpha= $' + str(np.round(np.abs(alpha[2]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[3], label=r'QS_t $\alpha= $' + str(np.round(np.abs(alpha[3]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#
##ax.plot(T, f_qs_r[0], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigmaD[0]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
##ax.plot(T, f_qs_r[1], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigmaD[1]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
##ax.plot(T, f_qs_r[2], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigmaD[2]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
##ax.plot(T, f_qs_r[3], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigmaD[3]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#
#
#clasical = np.ones_like(f_tmsv[0]) * .5
#plt.plot(T, clasical, 'r--')
#
#ax.set_xlabel(r'Fixed Loss')
#ax.set_ylabel(r'$\bar{\mathcal{F}} $')
## plt.ylim([-0.005, 0.2])
#
##ax.set_xlim(-0.1, 13)
#ax.set_ylim(.3, 1)
##ax.set_xscale('log')
#ax.grid()
#ax.legend()
#
#### To MATLAB
##data = [Tdb, f_tele1, f_tele2, f_tele3, f_tele4, f_dir1, f_dir2, f_dir3, f_dir4]
##scipy.io.savemat('matlab/fixed.mat', {'data':data})
#
#

############################# FIXED LOSS FIDELITIES PLOT
fig = plt.figure()
ax = fig.add_subplot(111)

T = np.linspace(0.01, 0.99, 15)
# T = [0.74, 0.81, 0.9, 1]
#T = np.logspace(-3, 0, 15)
eps = 0.005
eta = np.sqrt(10**(-1/10))
sigma = [2, 5, 10, 20]
# eta = 1

f_tmsv = []
f_qs_t = []
f_sb = []

for sig_ind in range(len(sigma)):
    f_tmsv_i = []
    f_qs_t_i = []
    f_sb_i = []

    for it in T:
        f_tmsv_i += [tmsv.opt_avg_fidelity(it, eps, eta, sigma[sig_ind])]
        f_qs_t_i += [qs_t.opt_avg_fidelity(it, eps, eta, sigma[sig_ind])]
#        f_qs_r_i += [qs_r.opt_avg_fidelity(it, eps, eta, sigmaT[sig_ind])]
        f_sb_i += [sb.opt_avg_fidelity(it, eps, np.sqrt(eta), sigma[sig_ind])]


    f_tmsv += [f_tmsv_i]
    f_qs_t += [f_qs_t_i]
    f_sb += [f_sb_i]




ax.plot(T, f_tmsv[0], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[0]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
ax.plot(T, f_tmsv[1], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[1]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
ax.plot(T, f_tmsv[2], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[2]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)
ax.plot(T, f_tmsv[3], label=r'TMSV $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='-', marker='o', markersize=4, linewidth=1.5)

#ax.plot(T, f_qs_t[0], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[0]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[1], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[1]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[2], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[2]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)
#ax.plot(T, f_qs_t[3], label=r'QS_t $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='-', marker='*', markersize=4, linewidth=1.5)

ax.plot(T, f_sb[0], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigma[0]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
ax.plot(T, f_sb[1], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigma[1]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
ax.plot(T, f_sb[2], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigma[2]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)
ax.plot(T, f_sb[3], label=r'QS_r $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='-', marker='v', markersize=4, linewidth=1.5)


clasical = np.ones_like(f_tmsv[0]) * .5
plt.plot(T, clasical, 'r--')

ax.set_xlabel(r'T')
ax.set_ylabel(r'$\bar{\mathcal{F}} $')
# plt.ylim([-0.005, 0.2])

#ax.set_xlim(-0.1, 13)
#ax.set_ylim(.49, 1)
#ax.set_xscale('log')
ax.grid()
ax.legend()
#ax.set_xscale('log')
### To MATLAB
data = [T, f_tmsv, f_qs_t, f_sb]
scipy.io.savemat('figs_noG/fixed_py.mat', {'data':data})

plt.show()
