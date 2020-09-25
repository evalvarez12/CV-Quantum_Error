# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 16:14:21 2020

@author: z5239621
"""

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


db = 0
eta = 10**(-0/10)
t_ext = 10**(-db/10)



colors_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors_fill = {'up':'darkblue', 'upEM':'skyblue', 'down':'darkgreen', 'downEM':'lawngreen'}
colors = {'uplink':'b', 'downlink':'g'}

#zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
zs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
rs = [1]

scint_up = {}
scint_down = {}

bins = np.arange(0.01, 1, .001)
colors = {0: colors_cycle[0], 15:colors_cycle[1], 30:colors_cycle[2], 45:colors_cycle[3], 60:colors_cycle[4], 75:colors_cycle[5]}

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
       data_file = "../../Laser_propagation/data/TELE_DOWN_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
       downlink_data = scipy.io.loadmat(data_file)['res'].transpose()

       data_file = "../../Laser_propagation/data/TELE_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_10000"
#       data_file = "../../Laser_propagation/data/TEST_UP_I_r=" + str(r) + "_z=" + str(z) + "_1024_50"
       uplink_data = scipy.io.loadmat(data_file)['res'].transpose()

       T_down = downlink_data * t_ext
       T_up = uplink_data * t_ext

       T_up = T_up.transpose()
       T_down = T_down.transpose()

       downlink_avg += [np.average(T_down)]
       downlink_std += [np.std(T_down)]

       uplink_avg += [np.average(T_up)]
       uplink_std += [np.std(T_up)]

       data_file = "../../Laser_propagation/data/DOWNLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       downlinkEM_data = scipy.io.loadmat(data_file)['t'].transpose()

       data_file = "../../Laser_propagation/data/UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
#       data_file = "../../Laser_propagation/data/TEST_UPLINK_PDF_EM_r=" + str(r) + '_z=' + str(z)
       uplinkEM_data = scipy.io.loadmat(data_file)['t'].transpose()

       T_downEM = downlinkEM_data * t_ext
       T_upEM = uplinkEM_data * t_ext

       downlinkEM_avg += [np.average(T_downEM)]
       downlinkEM_std += [np.std(T_downEM)]

       uplinkEM_avg += [np.average(T_upEM)]
       uplinkEM_std += [np.std(T_upEM)]

       print('z =', z, 'DOWN avg =', np.average(T_down))
       print('z =', z, 'UP avg =', np.average(T_up))

       print('---')

       print('z =', z, 'DOWN scint =', np.average(T_down**2)/(np.average(T_down)**2) - 1)
       print('z =', z, 'DOWN scint EM =', np.average(T_downEM**2)/(np.average(T_downEM)**2) - 1)

       print('z =', z, 'UP scint =', np.average(T_up**2)/(np.average(T_up)**2) - 1)
       print('z =', z, 'UP scint EM =', np.average(T_upEM**2)/(np.average(T_upEM)**2) - 1)


       data_file = "../../Laser_propagation/data/1550_UP_CENTER_z=" + str(z) + "_1024_10000"
#       data_file = "../../Laser_propagation/data/TEST_UP_CENTER_z=" + str(z) + "_1024_50"
       up_center_data = scipy.io.loadmat(data_file)['center'].transpose()


       data_file = "../../Laser_propagation/data/1550_DOWN_CENTER_z=" + str(z) + "_1024_10000"
#       data_file = "../../Laser_propagation/data/TEST_DOWN_CENTER_z=" + str(z) + "_1024_50"
       down_center_data = scipy.io.loadmat(data_file)['center'].transpose()

       T_up_center = up_center_data.transpose()
       T_down_center = down_center_data.transpose()

       print('---')

       print('z =', z, 'UP center scint =', np.average(T_up_center**2)/(np.average(T_up_center)**2) - 1)
       print('z =', z, 'DOWN center scint =', np.average(T_down_center**2)/(np.average(T_down_center)**2) - 1)

       print('---------------')


       scint_down.update( {z : np.average(T_down**2)/(np.average(T_down)**2) - 1} )
       scint_up.update({z : np.average(T_up**2)/(np.average(T_up)**2) - 1} )



       # plt.hist(T_upEM, histtype='step', label=str(z), density=True, linewidth=2, color=colors[z], linestyle=':')
       # plt.hist(T_up, histtype='step', label=str(z), density=True, linewidth=2, color=colors[z], linestyle='-')
