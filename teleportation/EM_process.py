# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:24:13 2019

@author: z5239621
"""

import TMSV as tmsv
import coherent as co
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"
############################# TRANSMISSIVITIES PLOT


# t_ext = 0.5012
db = 2
eta = 10**(-0/10)
t_ext = 10**(-db/10)




colors_fill = {'up':'darkblue', 'upEM':'skyblue', 'down':'darkgreen', 'downEM':'lawngreen'}
colors = {'uplink':'b', 'downlink':'g'}

zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
rs = [0.75]
fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax2 = ax1.twinx()

scint_up = {}
scint_down = {}

for r in rs:
   uplink_avg = []
   uplink_std = []
   downlink_avg = []
   downlink_std = []
   uplinkEM_avg = []
   uplinkEM_std = []
   downlinkEM_avg = []
   downlinkEM_std = []

   for z in zs:

       tag = "_z=" + str(z) + "_1024_10000"
       data_file = "data/1550_DOWN_I_r=" + str(r) + tag
       downlink_data = scipy.io.loadmat(data_file)['res'].transpose()

       data_file = "data/1550_UP_I_r=" + str(r) + tag
       uplink_data = scipy.io.loadmat(data_file)['res'].transpose()

       T_down = downlink_data * t_ext
       T_up = uplink_data * t_ext

       downlink_avg += [np.average(T_down)]
       downlink_std += [np.std(T_down)]

       uplink_avg += [np.average(T_up)]
       uplink_std += [np.std(T_up)]

       data_file = "data/DOWNLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       downlinkEM_data = scipy.io.loadmat(data_file)['t'].transpose()

       data_file = "data/UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       uplinkEM_data = scipy.io.loadmat(data_file)['t'].transpose()

       T_downEM = downlinkEM_data * t_ext
       T_upEM = uplinkEM_data * t_ext

       downlinkEM_avg += [np.average(T_downEM)]
       downlinkEM_std += [np.std(T_downEM)]

       uplinkEM_avg += [np.average(T_upEM)]
       uplinkEM_std += [np.std(T_upEM)]


       print('z =', z, 'DOWN scint =', np.average(T_down**2)/(np.average(T_down)**2) - 1)
       print('z =', z, 'DOWN scint EM =', np.average(T_downEM**2)/(np.average(T_downEM)**2) - 1)

       print('z =', z, 'UP scint =', np.average(T_up**2)/(np.average(T_up)**2) - 1)
       print('z =', z, 'UP scint EM =', np.average(T_upEM**2)/(np.average(T_upEM)**2) - 1)
       print('---------------')


       scint_down.update( {z : np.average(T_down**2)/(np.average(T_down)**2) - 1} )
       scint_up.update({z : np.average(T_up**2)/(np.average(T_up)**2) - 1} )


   # Y1 = [y - e for y, e in zip(Tavg, Tstd)]
   # Y2 = [y + e for y, e in zip(Tavg, Tstd)]
   #
   # Y1_EM = [y - e for y, e in zip(Tavg_EM, Tstd_EM)]
   # Y2_EM = [y + e for y, e in zip(Tavg_EM, Tstd_EM)]
   #
   # ax1.fill_between(zs, Y1_EM, Y2_EM, alpha=0.25, color=colors_fill[str(r)+'EM'])
   # ax1.errorbar(zs, Tavg_EM, Tstd_EM, label=r'$r_d=' + str(r)+'$m EM', alpha=1, fmt=':', capsize=4, markersize=4, linewidth=1.5, c=colors_fill[str(r)+'EM'])
   #
   # ax1.fill_between(zs, Y1, Y2, alpha=0.5, color=colors_fill[str(r)])
   # ax1.errorbar(zs, Tavg, Tstd, label=r'$r_d=' + str(r)+'$m', alpha=1, fmt=':', capsize=4, markersize=4, linewidth=1.5, c=colors_fill[str(r)])


   ax1.errorbar(zs, downlink_avg, downlink_std, label=r'downlink', linestyle='-', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['downlink'])
   ax1.errorbar(zs, uplink_avg, uplink_std, label=r'uplink', linestyle='-', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['uplink'])

   ax1.errorbar(zs, downlinkEM_avg, downlinkEM_std, label=r'downlink EM', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['downlink'])
   ax1.errorbar(zs, uplinkEM_avg, uplinkEM_std, label=r'uplink EM', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['uplink'])



ax1.set_xlabel(r'$\zeta$ (deg)')
ax1.set_ylabel(r'$T$')
# ax2.set_ylabel(r'$\gamma$')
ax1.grid()
ax1.legend()









plt.figure()
#for z in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]:
r = 0.75

colors_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = {0: colors_cycle[0], 15:colors_cycle[1], 30:colors_cycle[2], 45:colors_cycle[3], 60:colors_cycle[4], 75:colors_cycle[5]}

bins = np.arange(0.03, .450, .001)
bins2 = np.arange(0.0, 1, .001)

for z in [0, 15, 30, 45, 60]:
    data_file = "data/DOWNLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
    downlink_data = scipy.io.loadmat(data_file)['t'].transpose()

    data_file = "data/UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
    uplink_data = scipy.io.loadmat(data_file)['t'].transpose()

    print('z = ', z)
    plt.hist(uplink_data * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=2, color=colors[z], linestyle=':')

#    plt.hist((1-P)/P * db2, bins=bins2, histtype='step', label=str(z), density=True, linewidth=1.5, color=colors[z])


# for z in [0, 30, 75]:
#     tag = "_z=" + str(z) + "_1024_10000"
#     data_file = "l0d=0.1SH_DOWN_I_r=" + str(r)  + tag
#     I = scipy.io.loadmat(data_file)['res']
#     plt.hist(I * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=1.5, linestyle=':', color=colors[z])


#    plt.hist(P0*I, bins=bins)
    #plt.hist(I, bins=bins)
plt.xlabel('$T$')
plt.ylabel('PDF')
plt.legend()


colors = {'uplink':'b', 'downlink':'g'}
zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
rs = [0.75]
fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax2 = ax1.twinx()
for r in rs:
   uplink_avg = []
   uplink_std = []
   downlink_avg = []
   downlink_std = []
   for z in zs:

       data_file = "data/DOWNLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       downlink_data = scipy.io.loadmat(data_file)['t'].transpose()

       data_file = "data/UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       uplink_data = scipy.io.loadmat(data_file)['t'].transpose()


       T_down = downlink_data * t_ext
       T_up = uplink_data * t_ext

       downlink_avg += [np.average(T_down)]
       downlink_std += [np.std(T_down)]

       uplink_avg += [np.average(T_up)]
       uplink_std += [np.std(T_up)]

   # ax2.errorbar(zs, Pavg, Pstd,  linestyle=':', label=r'$r_d=' + str(r)+'$m', marker='^', capsize=1, markersize=4, c=colorsP[r])

   ax1.errorbar(zs, downlink_avg, downlink_std, label=r'downlink', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['downlink'])
   ax1.errorbar(zs, uplink_avg, uplink_std, label=r'uplink', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c=colors['uplink'])


   ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

   color = 'tab:red'
   ax2.set_ylabel(r'$\sigma_{SI}$', color=color)  # we already handled the x-label with ax1
   ax2.plot(zs, list(scint_down.values())[:-1], 'v-', label=r'downlink', color='red')
   ax2.plot(zs, list(scint_up.values())[:-1], 'v-', label=r'uplink', color='orange')
   ax2.tick_params(axis='y', labelcolor=color)
   ax2.legend()


ax1.set_xlabel(r'$\zeta$ (deg)')
ax1.set_ylabel(r'$T$')
# ax2.set_ylabel(r'$\gamma$')
ax1.grid()

ax1.legend()


############################## FIDELITIES PLOT
fig = plt.figure()
ax1 = fig.add_subplot(111)

zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
rs = [0.75]

eps1 = 0.02
eps2 = 0.07

alpha = [(1 + 1j)*.8, 1 + 1j, (1 + 1j)*1.2,]


for r in rs:
   uplink_avg = []
   uplink_std = []
   downlink_avg = []
   downlink_std = []

   f_up1_avg = []
   f_up1_std = []
   f_up2_avg = []
   f_up2_std = []
   f_up3_avg = []
   f_up3_std = []

   f_tmsv_avg = []
   f_tmsv_std = []

   for z in zs:

       data_file = "data/DOWNLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       downlink_data = scipy.io.loadmat(data_file)['t'].transpose()

       data_file = "data/UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       uplink_data = scipy.io.loadmat(data_file)['t'].transpose()


       T_down = downlink_data * t_ext
       T_up = uplink_data * t_ext

       downlink_avg += [np.average(T_down)]
       downlink_std += [np.std(T_down)]

       uplink_avg += [np.average(T_up)]
       uplink_std += [np.std(T_up)]

       # TMSV
       V_opt, F_opt = tmsv.opt_values(np.average(T_down), eps1, eta, alpha=1 + 1j)
       print("z =", z, V_opt, F_opt)

       eps = scint_down[z] * V_opt
       f_tmsv_avg += [np.average(tmsv.fidelity(V_opt, T_down, eps, eta, alpha))]
       f_tmsv_std += [np.std(tmsv.fidelity(V_opt, T_down, eps, eta, alpha))]

       eps = scint_up[z] * np.abs(alpha[0])
       f_up1_avg += [np.average(co.fidelity(T_up, eps, alpha[0]))]
       f_up1_std += [np.std(co.fidelity(T_up, eps, alpha[0]))]

       eps = scint_up[z] * np.abs(alpha[1])
       f_up2_avg += [np.average(co.fidelity(T_up, eps, alpha[1]))]
       f_up2_std += [np.std(co.fidelity(T_up, eps, alpha[1]))]

       eps = scint_up[z] * np.abs(alpha[2])
       f_up3_avg += [np.average(co.fidelity(T_up, eps, alpha[2]))]
       f_up3_std += [np.std(co.fidelity(T_up, eps, alpha[2]))]

   ax1.errorbar(zs, f_tmsv_avg, f_tmsv_std, label=r'Downlink', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
   ax1.errorbar(zs, f_up1_avg, f_up1_std, label=r'Uplink $|\alpha|^2= $' + str(np.round(np.abs(alpha[0]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
   ax1.errorbar(zs, f_up2_avg, f_up2_std, label=r'Uplink $|\alpha|^2= $' + str(np.round(np.abs(alpha[1]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
   ax1.errorbar(zs, f_up3_avg, f_up3_std, label=r'Uplink $|\alpha|^2= $' + str(np.round(np.abs(alpha[2]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)

#plt.plot(zs, np.zeros_like(zs), ls='--', c='k')

clasical = np.ones_like(f_tmsv_avg) * .5
plt.plot(zs, clasical, 'r--')


ax1.set_xlabel(r'$\zeta$ (deg)')
ax1.set_ylabel(r'Fidelity')
# plt.ylim([-0.005, 0.2])
ax1.grid()
ax1.legend()





plt.show()
