# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:38:54 2020

@author: z5239621
"""

import scipy.io
import numpy as np
import matplotlib.pyplot as plt


############ UPLINK & DOWNLINK CHANNELS PROPERTIES PLOT

db = 0
eta = 10**(-0/10)
t_ext = 10**(-db/10)


scint_up = {}
scint_down = {}

colors = {'uplink':'b', 'downlink':'g'}
zs = [0, 10, 15, 20 , 25]
rs = [1]
fig = plt.figure()
ax = fig.add_subplot(111)
# ax2 = ax.twinx()
for r in rs:
   uplink_avg = []
   uplink_std = []
   downlink_avg = []
   downlink_std = []
   for z in zs:

       data_file = "../../Laser_propagation/data/TELE_DOWN_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
       downlink_data = scipy.io.loadmat(data_file)['res'].transpose()

#       data_file = "../../Laser_propagation/data/1550_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
       data_file = "../../Laser_propagation/data/TELE_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
       uplink_data = scipy.io.loadmat(data_file)['res'].transpose()


       T_down = downlink_data * t_ext
       T_up = uplink_data * t_ext
       
#       T_down = -10 * np.log10(T_down)
#       T_up = -10 * np.log10(T_up)

       downlink_avg += [np.average(T_down)]
       downlink_std += [np.std(T_down)]

       uplink_avg += [np.average(T_up)]
       uplink_std += [np.std(T_up)]
       
       scint_down.update( {z : np.average(T_down**2)/(np.average(T_down)**2) - 1} )
       scint_up.update({z : np.average(T_up**2)/(np.average(T_up)**2) - 1} )

   # ax2.errorbar(zs, Pavg, Pstd,  linestyle=':', label=r'$r_d=' + str(r)+'$m', marker='^', capsize=1, markersize=4, c=colorsP[r])

   ax.errorbar(zs, downlink_avg, downlink_std, label=r'downlink', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c='navy')
   ax.errorbar(zs, uplink_avg, uplink_std, label=r'uplink', linestyle='--', marker='s', capsize=4, markersize=4, linewidth=1.5, c='dodgerblue')


   ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

   ax2.set_ylabel(r'$\sigma_{SI}$', color='red')  # we already handled the x-label with ax
   ax2.plot(zs, list(scint_down.values()), 'v-', label=r'downlink', color='crimson')
   ax2.plot(zs, list(scint_up.values()), '^-', label=r'uplink', color='tomato')
   ax2.tick_params(axis='y', labelcolor='red')
   ax2.set_ylim(0.035, -0.001)
   ax2.legend()

#ax.set_zorder(1)  
#ax.patch.set_visible(False)
ax.set_xlabel(r'$\zeta$ (deg)')
ax.set_ylabel(r'$\langle t \rangle (dB) $', color='blue')
ax.tick_params(axis='y', labelcolor='blue')
# ax2.set_ylabel(r'$\gamma$')
ax.grid()
#ax.set_ylim(2.9, 17.5)

ax.legend()

plt.show()