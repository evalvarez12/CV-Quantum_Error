# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:38:54 2020

@author: z5239621
"""

import scipy.io
import numpy as np
import matplotlib.pyplot as plt


############ UPLINK & DOWNLINK CHANNELS PROPERTIES PLOT

colors = {'uplink':'b', 'downlink':'g'}
fig = plt.figure()
ax = fig.add_subplot(111)
uplink_avg = []
uplink_std = []
downlink_avg = []
downlink_std = []
r = 1
bins = np.arange(0., .55, .01)

for z in [10]:

   data_file = "../../Laser_propagation/data/TELE__L0inf_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
   downlink_data = scipy.io.loadmat(data_file)['res'].transpose()[0]

#       data_file = "../../Laser_propagation/data/1550_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
   data_file = "../../Laser_propagation/data/TELE_L0inf2_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
   uplink_data = scipy.io.loadmat(data_file)['res'].transpose()[0]

   downlink_avg += [np.average(downlink_data)]
   downlink_std += [np.std(downlink_data)]

   uplink_avg += [np.average(uplink_data)]
   uplink_std += [np.std(uplink_data)]


   ax.hist(uplink_data, bins=bins, histtype='step', label='uplink z-'+str(z), density=True, linewidth=2, linestyle='-')
   ax.hist(downlink_data, bins=bins, histtype='step', label='downlink z-'+str(z), density=True, linewidth=2, linestyle='-')

# ax.errorbar(zs, downlink_avg, downlink_std, label=r'downlink', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5, c='navy')
# ax.errorbar(zs, uplink_avg, uplink_std, label=r'uplink', linestyle='--', marker='s', capsize=4, markersize=4, linewidth=1.5, c='dodgerblue')


# ax.set_xlabel(r'$\zeta$ (deg)')
# ax.set_ylabel(r'$\langle t \rangle$', color='blue')
# ax.tick_params(axis='y', labelcolor='blue')
ax.grid()
ax.legend()
ax.set_title(r'$\omega_0=15cm$')
plt.show()

print('AVG DOWN:', downlink_avg)
print('AVG UP:', uplink_avg)
print('STD DOWN:', downlink_std)
print('STD UP:', uplink_std)

