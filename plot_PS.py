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


def plot(N, params, MG=True):
    key_rates = []
    indeces = []
    options = ['none', 'tps', 'rps', 'tsc', 'rsc', 'tct', 'rct']
#    options = ['none', 'tps', 'rps']

    for option in options:

        if option == 'rsc':
            filename = names.measurements_line(10, 'KR', params, option)
            filename_ind = names.indeces_line(10, 'KR', params, option, 'eta')
        else:
            filename = names.measurements_line(N, 'KR', params, option)
            filename_ind = names.indeces_line(N, 'KR', params, option, 'eta')


        k_rate = np.load(filename + ".npy")
        key_rates += [k_rate]

        inds = np.load(filename_ind + '.npy')
        indeces += [inds]


    fig1, ax = plt.subplots()
    lines_types = ['k*-', 'b*-', 'r*-', 'y*-', 'm*-', 'c*-', 'g*-']
    for i in range(len(key_rates)):
        ax.plot(indeces[i], key_rates[i], lines_types[i], label=options[i])


    if MG:
        mg_ratesNO = scipy.io.loadmat('data/rate_no.mat')
        mg_ratesR = scipy.io.loadmat('data/rate_r.mat')
        mg_ratesT = scipy.io.loadmat('data/rate_t.mat')
        mg_indNO = scipy.io.loadmat('data/ind_no.mat')
        mg_indR = scipy.io.loadmat('data/ind_r.mat')
        mg_indT = scipy.io.loadmat('data/ind_t.mat')

        mg_ratesNO = mg_ratesNO['RateNo'][0]
        mg_ratesR = mg_ratesR['RateR'][0]
        mg_ratesT = mg_ratesT['RateT'][0]
        mg_indNO = mg_indNO['indd'][0]
        mg_indR = mg_indR['indR'][0]
        mg_indT = mg_indT['indT'][0]

        ax.plot(mg_indNO, mg_ratesNO, 'k--')
        ax.plot(mg_indR, mg_ratesR, 'r--')
        ax.plot(mg_indT, mg_ratesT, 'b--')

    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("Key rate")
    #ax.legend(["No PS", "Transmitter PS", "Receiver PS", "Transmitter SC", "Mg no PS", "Mg r-PS", "Mg t-PS"])


    #key_rates1b = np.load(filename1 + "BAD.npy")
    #key_rates2b = np.load(filename2 + "BAD.npy")
    #key_rates3b = np.load(filename3 + "BAD.npy")
    #ax.plot(tes, key_rates1b, 'ko')
    #ax.plot(tes, key_rates2b, 'bo')
    #ax.plot(tes, key_rates3b, 'ro')

    ax.set_yscale('log')


    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
    ax.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    #ax.minorticks_on()
    plt.legend()
    plt.show()



def plot_split(N, params, MG=True):
    i_stolens = []
    i_shareds = []
    indeces = []
    options = ['none', 'tps', 'rps', 'tsc', 'rsc']
#    options = ['none', 'tps', 'rps']

    for option in options:

        if option == 'rsc':
            filename1 = names.measurements_line(10, 'I_shared', params, option)
            filename2 = names.measurements_line(10, 'I_stolen', params, option)
            filename_ind = names.indeces_line(10, 'I_s', params, option, 'eta')
        else:
            filename1 = names.measurements_line(N, 'I_shared', params, option)
            filename2 = names.measurements_line(N, 'I_stolen', params, option)
            filename_ind = names.indeces_line(N, 'I_s', params, option, 'eta')


        i_shared = np.load(filename1 + ".npy")
        i_shareds += [i_shared]

        i_stolen = np.load(filename2 + ".npy")
        i_stolens += [i_stolen]

        inds = np.load(filename_ind + '.npy')
        indeces += [inds]


    fig1, ax = plt.subplots()
    lines_types = ['k*-', 'b*-', 'r*-', 'y*-', 'm*-', 'c*-', 'g*-']
    lines_types2 = ['ko-', 'bo-', 'ro-', 'yo-', 'mo-', 'co-', 'go-']

    for i in range(len(i_shareds)):
        ax.plot(indeces[i], i_shareds[i], lines_types[i], label=options[i])
        ax.plot(indeces[i], i_stolens[i], lines_types2[i], label=options[i])

    print(params)
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("Key rate")
    #ax.legend(["No PS", "Transmitter PS", "Receiver PS", "Transmitter SC", "Mg no PS", "Mg r-PS", "Mg t-PS"])


    #key_rates1b = np.load(filename1 + "BAD.npy")
    #key_rates2b = np.load(filename2 + "BAD.npy")
    #key_rates3b = np.load(filename3 + "BAD.npy")
    #ax.plot(tes, key_rates1b, 'ko')
    #ax.plot(tes, key_rates2b, 'bo')
    #ax.plot(tes, key_rates3b, 'ro')

    ax.set_yscale('log')


    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
    ax.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    #ax.minorticks_on()
    plt.legend()
    plt.show()





## Parameters
#N = 20
#mpn = 1.3
#mpne = 0.001
#f = 0.95
#option = 'tsc'
## Operations options
#t = .1
#k_ps = 0.9
#k_sc = 0.1
#k_ct = 0.5
#
##params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "t=" + str(t)]
#params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "k_ps=" + str(k_ps),
#          "k_sc=" + str(k_sc), "k_ct=" + str(k_ct)]
#
#plot(N, params)
#
#
# Parameters
N = 20
mpn = 1.3
mpne = 0.001
f = 0.95
option = 'tsc'
# Operations options
t = .1
k_ps = 0.9
k_sc = 0.1
k_ct = 0.5

#params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "t=" + str(t)]
params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "k_ps=" + str(k_ps),
         "k_sc=" + str(k_sc), "k_ct=" + str(k_ct)]

plot_split(N, params)
