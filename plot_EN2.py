# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:04:09 2019

@author: Eduardo Villasenor
"""

import src.names as names
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.io


def plot(N, params):
    datas = []
    datas2 = []
    indeces = []
    options = ['none', 'tps', 'rps', 'tsc', 'rsc']
#    options = ['none', 'tps', 'rps']

    for option in options:

        filename = names.measurements_line(N, 'EN2', params, option)
        filename = names.measurements_line(N, 'EN2_p', params, option)
        filename_ind = names.indeces_line(N, 'EN2', params, option, 'eta')

        data = np.load(filename + ".npy")
        datas += [data]
        
        data2 = np.load(filename + ".npy")
        datas2 += [data2]

        inds = np.load(filename_ind + '.npy')
        indeces += [inds]


    fig1, ax = plt.subplots()
    lines_types = ['k*-', 'b*-', 'r*-', 'y*-', 'm*-', 'c*-', 'g*-']
    for i in range(len(indeces)):
        
        x = -10*np.log10(indeces[i])
        
#        y = np.multiply(datas[i], datas2[i])
        y = datas[i]
        
        ax.plot(x, y, lines_types[i], label=options[i])


    ax.set_xlabel(r"Atenuation (dB)")
    ax.set_ylabel("$E_N$")
    #ax.legend(["No PS", "Transmitter PS", "Receiver PS", "Transmitter SC", "Mg no PS", "Mg r-PS", "Mg t-PS"])


    #key_rates1b = np.load(filename1 + "BAD.npy")
    #key_rates2b = np.load(filename2 + "BAD.npy")
    #key_rates3b = np.load(filename3 + "BAD.npy")
    #ax.plot(tes, key_rates1b, 'ko')
    #ax.plot(tes, key_rates2b, 'bo')
    #ax.plot(tes, key_rates3b, 'ro')

#    ax.set_yscale('log')


#    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
#    ax.yaxis.set_major_locator(locmaj)
#    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
#    ax.yaxis.set_minor_locator(locmin)
#    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    #ax.minorticks_on()
    plt.legend()
    plt.show()





## Parameters
#N = 5
#mpn = 0.1
## Operations options
#k_ps = 0.5
#k_sc = 0.01
#k_ct = 0.5
#
##params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "t=" + str(t)]
#params = ["mpn=" + str(mpn), "k_ps=" + str(k_ps),
#          "k_sc=" + str(k_sc), "k_ct=" + str(k_ct)]
#
#plot(N, params)
#
#
