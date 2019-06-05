# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:04:09 2019

@author: Eduardo Villasenor
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.io

def plot():
    key_rates = []
    indeces = []
    for ext in ['none', 'tps', 'rps', 'tsc_e', 'rsc_e', 'tct', 'rct']:
        f_name = "data/result_PS_" + ext
        k_rate = np.load(f_name + ".npy")
        key_rates += [k_rate]
    
        i_name = "data/indeces_PS_" + ext
        inds = np.load(i_name + '.npy')
        indeces += [inds]
    
    
    fig1, ax = plt.subplots()
    lines_types = ['k*-', 'b*-', 'r*-', 'y*-', 'm*-', 'c*-', 'g*-']
    for i in range(len(key_rates)):
        ax.plot(indeces[i], key_rates[i], lines_types[i])
    
    
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
    plt.show()
