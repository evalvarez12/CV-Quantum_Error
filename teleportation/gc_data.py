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
#import matplotlib.pyplot as plt



num = 30
Ts = np.linspace(0.0001, 1, num)
es = np.logspace(-2.2, -0.9, num)

eta = np.sqrt(10**(-1/10))
sigma = 10
# eta = 1
#eta = 0.7

f_tmsv = np.zeros([num, num])
f_sb = np.zeros([num, num])

for i in range(num):
    for j in range(num):
        T = Ts[i]
        e = es[j]
        
        f_tmsv[i, j] = tmsv.opt_avg_fidelity(T, e, eta, sigma)
        f_sb[i, j] = sb.opt_avg_fidelity(T, e, eta, sigma)


### To MATLAB
data = [f_tmsv, f_sb]
scipy.io.savemat('figs_noG/data-py-sig10.mat', {'data':data})

