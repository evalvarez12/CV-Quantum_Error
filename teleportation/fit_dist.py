# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:38:02 2020

@author: z5239621
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import scipy.io


t_ext = 0.5012

plt.figure()
#for z in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]:
r = 0.75
si = []
for z in [0, 5 ,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]:

    tag = "_z=" + str(z) + "_1024_10000"


    data_file = "DOWN_I_r=" + str(r)  + tag
    I = scipy.io.loadmat(data_file)['res']
    data_file = "DOWN_P0_r=" + str(r)  + tag
    P0 = scipy.io.loadmat(data_file)['snr0']
    data_file = "DOWN_P_r=" + str(r)  + tag
    P = scipy.io.loadmat(data_file)['snr']


    T = I * t_ext

    print("z=", z, len(P))
    print("var=", np.var(T))
    print("SI=", np.var(T)/(np.mean(T)**2)* 10)
    bins = np.arange(0.1, .5, .0005)
    plt.hist(I * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=1.5)
#    plt.hist(P0*I, bins=bins)
    #plt.hist(I, bins=bins)
    
    si = si + [np.var(T)/(np.mean(T)**2)]
    
plt.xlabel('$T$')
plt.ylabel('PDF')
#plt.legend()


#plt.figure()
#
#zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
#
#plt.plot(zs, si, 'bo')
#plt.xlabel('$\zeta$')
#plt.ylabel(r'$\tilde{\sigma}(T)$')
#plt.yscale('log')
#
#plt.figure()
#
#
#colors = {1:'b', 0.75:'g'}
#
#zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
#rs = [0.75, 1]
#for r in rs:
#    Tavg = []
#    Tstd = []
#    for z in zs:
#
#        tag = "_z=" + str(z) + "_1024_10000"
#        data_file = "DOWN_I_r=" + str(r) + tag
#        I = scipy.io.loadmat(data_file)['res']
#        data_file = "DOWN_P0_r=" + str(r) + tag
#        P0 = scipy.io.loadmat(data_file)['snr0']
#        data_file = "DOWN_P_r=" + str(r) + tag
#        P = scipy.io.loadmat(data_file)['snr']
#
#        I = I.transpose()[0]
#        P = P.transpose()[0]
#        P0 = P0.transpose()[0]
#
#        T =  I * t_ext
#
#
#        Tavg += [np.average(T)]
#        Tstd += [np.std(T)]
#
#    plt.errorbar(zs, Tavg, Tstd, label=r'$r_d=' + str(r)+'$m', marker='s', capsize=4, markersize=4, c=colors[r])
#plt.xlabel(r'$\zeta$ (deg)')
#plt.ylabel(r'$T$')
#plt.grid()
#
#plt.legend()


########################################################  Compare with other 10d data with l0 = 0.005 * L0
#r =0.75
#z = 30
#data_file1 = "l0d_DOWN_I_r=" + str(r)  + "_z=" + str(z) + "_1024_10000"
#I1 = scipy.io.loadmat(data_file1)['res']
#
#plt.hist(I1 * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=3)
#
#z = 75
#data_file2 = "l0d_DOWN_I_r=" + str(r)  + "_z=" + str(z) + "_1024_10000"
#I2 = scipy.io.loadmat(data_file2)['res']
#plt.hist(I2 * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=3)
#############################################################################
