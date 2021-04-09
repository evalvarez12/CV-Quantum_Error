# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:43:48 2021

@author: z5239621
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

mat = scipy.io.loadmat('figs_noG/data-py-sig10.mat')

tmsv = mat['data'][0]
sb = mat['data'][1]

num = 30
Ts = np.linspace(0.0001, 1, num)
es = np.logspace(-2.2, -0.9, num)

plt.contour(es, Ts, sb)
plt.contour(es, Ts, tmsv, colors='k')

plt.xscale('log')
plt.colorbar()


plt.show()