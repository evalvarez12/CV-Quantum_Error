# -*- coding: utf-8 -*-
"""
Create results for the GLOBECOMM paper
Created on Mon Mar  8 11:32:18 2021

@author: z5239621
"""


import TMSV as tmsv
import TMSV_sq as tmsvsq
import qsc as qs_t
import qsc_rec as qs_r
import sb as sb
import scipy.io
import numpy as np
#import matplotlib.pyplot as plt



#num = 30
#Ts = np.linspace(0.0001, 1, num)
#es = np.logspace(-2.2, -0.9, num)
#
#eta = np.sqrt(10**(-1/10))
#sigma = 10
## eta = 1
##eta = 0.7
#
#f_tmsv = np.zeros([num, num])
#f_sb = np.zeros([num, num])
#
#for i in range(num):
#    for j in range(num):
#        T = Ts[i]
#        e = es[j]
#
#        f_tmsv[i, j] = tmsv.opt_avg_fidelity(T, e, eta, sigma)
#        f_sb[i, j] = sb.opt_avg_fidelity(T, e, eta, sigma)
#
#
#### To MATLAB
#data = [f_tmsv, f_sb]
#scipy.io.savemat('figs_noG/data-py-sig10.mat', {'data':data})




# Fiber
Ls = np.linspace(30, 150, 24);
num = len(Ls)
me = (0.0086-.0033)/100


# Sat
#Ts = [0.06143690562983698, 0.06086526299427021, 0.05917661503498977, 0.056443389085287135, 0.05323006348104734, 0.049082101315961746, 0.044150186248965786,  0.037171035021636226, 0.031114584577009226, 0.025515914750767258, 0.0198735248725103, 0.013815858165353754, 0.009422950626114226, 0.00582421, 0.00264689];
#num = 15

eta = np.sqrt(10**(-1/10))
#sq
# sigma = 1
# tmsv
sigma = 10
# eta = 1
#eta = 0.7

f_tmsv = np.zeros(num)



for i in range(num):

        # Fiber
        L = Ls[i]
        T = 10**(-(0.16 * L)/10);
        e = me * L + 0.0006

        # Sat
#        T = Ts[i];
#        e = T * 0.0186 + 0.0133/0.95;

        f_tmsv[i] = tmsv.opt_avg_fidelity(T, e, eta, sigma)
        # f_tmsv[i] = tmsvsq.opt_avg_fidelity(T, e, eta, sigma)

#        f_sb[i, j] = sb.opt_avg_fidelity(T, e, eta, sigma)




### To MATLAB
data = f_tmsv
scipy.io.savemat('figs_noG/py_fiber_sig10-2.mat', {'data':data})
