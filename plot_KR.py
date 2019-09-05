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


def plotMG(N, params, MG=True):
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
    
    
    
def plot(N, params):
    key_rates = []
    indeces = []
#    options = ['none', 'tps', 'rps', 'tqs', 'rqs']
    options = ['none', 'tps', 'rps', 'tqs', 'rqs', 'tot', 'rot', 'tct', 'rct', 'tpa', 'rpa']
#    options = ['none', 'tqs', 'rqs']

    for option in options:

        if option == 'rqs':
            filename = names.measurements_line(10, 'KR2', params, option)
            filename_ind = names.indeces_line(10, 'KR2', params, option, 'eta')
        else:
            filename = names.measurements_line(N, 'KR2', params, option)
            filename_ind = names.indeces_line(N, 'KR2', params, option, 'eta')


        k_rate = np.load(filename + ".npy")
        key_rates += [k_rate]

        inds = np.load(filename_ind + '.npy')
        indeces += [inds]


    fig1, ax = plt.subplots()
    
    ax.yaxis.grid(color='gray', linestyle='solid')
    ax.xaxis.grid(color='gray', linestyle='solid')
    
    lines_types = ['k*-', 'b*-', 'r*-', 'g*-', 'm*-', 'c*-', 'y*-', 'bs-', 'rs-', 'ys-', 'ms-']
    legends = ['TMSV', 'tps', 'rps', 'tqs', 'rqs', 'tot', 'rot', 'tct', 'rct', 'tpa', 'rpa']
    for i in range(len(key_rates)):
        x = -10*np.log10(indeces[i])
        ax.plot(x, key_rates[i], lines_types[i], label=legends[i], linewidth=3)


    plt.xlim(-0.1, 22.5)
#    plt.xlim(-0.1, 13.2)

    plt.rcParams["font.family"] = "Times New Roman"

    plt.title(r"Squeezing = 8 dB", size=15)
    ax.set_xlabel(r"Attenuation (dB)", size=15)
    ax.set_ylabel("Key rate (bit/pulse)", size=15)
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
    ax.tick_params(axis='both', which='major', labelsize=13)
    #ax.minorticks_on()
#    plt.legend(loc='upper right', fontsize=15)
    plt.show()



def plotIaltern(N, params):
    key_rates = []
    indeces = []
    
    key_ratesIaltern = []
    indecesIaltern = []
    
    options = ['none', 'tps', 'rps', 'tqs', 'rqs']
#    options = ['none', 'tps', 'rps', 'tqs', 'rqs', 'tot', 'rot', 'tct', 'rct', 'tpa', 'rpa']
#    options = ['none', 'tqs', 'rqs']

    for option in options:

        if option == 'rqs':
            filename = names.measurements_line(10, 'KR2', params, option)
            filename_ind = names.indeces_line(10, 'KR2', params, option, 'eta')
            
            filenameIaltern = names.measurements_line(10, 'KRIaltern', params, option)
            filename_indIaltern = names.indeces_line(10, 'KRIaltern', params, option, 'eta')
        else:
            filename = names.measurements_line(N, 'KR2', params, option)
            filename_ind = names.indeces_line(N, 'KR2', params, option, 'eta')
            
            filenameIaltern = names.measurements_line(N, 'KRIaltern', params, option)
            filename_indIaltern = names.indeces_line(N, 'KRIaltern', params, option, 'eta')


        k_rate = np.load(filename + ".npy")
        key_rates += [k_rate]

        inds = np.load(filename_ind + '.npy')
        indeces += [inds]


        k_rateIaltern = np.load(filenameIaltern + ".npy")
        key_ratesIaltern += [k_rateIaltern]

        indsIaltern = np.load(filename_indIaltern + '.npy')
        indecesIaltern += [indsIaltern]

    fig1, ax = plt.subplots()
    
    ax.yaxis.grid(color='gray', linestyle='solid')
    ax.xaxis.grid(color='gray', linestyle='solid')
    
    lines_types = ['k*-', 'b*-', 'r*-', 'g*-', 'm*-', 'c*-', 'y*-', 'bs-', 'rs-', 'ys-', 'ms-']
    legends = ['TMSV', 'tps', 'rps', 'tqs', 'rqs', 'tot', 'rot', 'tct', 'rct', 'tpa', 'rpa']
    for i in range(len(key_rates)):
        x = -10*np.log10(indeces[i])
        ax.plot(x, key_rates[i], lines_types[i], label=legends[i], linewidth=3)

        x = -10*np.log10(indecesIaltern[i])
        ax.plot(x, key_ratesIaltern[i], lines_types[i][0] + 's-', label=legends[i], linewidth=3)

    plt.xlim(-0.1, 22.5)
#    plt.xlim(-0.1, 13.2)

    plt.rcParams["font.family"] = "Times New Roman"

    plt.title(r"Squeezing = 8 dB", size=15)
    ax.set_xlabel(r"Attenuation (dB)", size=15)
    ax.set_ylabel("Key rate (bit/pulse)", size=15)
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
    ax.tick_params(axis='both', which='major', labelsize=13)
    #ax.minorticks_on()
#    plt.legend(loc='upper right', fontsize=15)
    plt.show()


######## EXAMPLE 
# Parameters
N = 20
r = .92
r_eve = 0.033

# Operations options
k_ps = 0.95
k_qs = 0.05
params = ["r=" + str(r), "r_eve=" + str(r_eve), "k_ps=" + str(k_ps),
          "k_qs=" + str(k_qs)]

plot(N, params)