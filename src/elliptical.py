# -*- coding: utf-8 -*-
"""
Object to simulate the atmospheric effects over CV states
Created on Thu Sep  5 14:37:49 2019

@author: Eduardo Villasenor
"""

import numpy as np

def theta(sci, sigma, size):
    a  = sci**2 * sigma**(5/6)
    mean = np.log((1+2.96 * a)**2 / (sigma**2 * np.sqrt((1 + 2.96 * a)**2 + 1.2 * a)))
    
    var = np.log(1 + (1.2 * a) / (1 + 2.96 * a)**2)
#    var = var - mean**2
    
    var12 = np.log(1 - (0.8 * a) / (1 + 2.96 * a)**2)
#    var12 = var12 - mean**2
    
    mean_vec = [mean, mean]
    var_vec = np.array([[var, var12], [var12, var]])
    print(mean)
    print( var_vec)
    return np.random.multivariate_normal(mean_vec, var_vec, size)
    
    
    

k = 1/(1550 * 1e-7)
w0 = 80
a = 110
L = 10 * 1e6

sigma = (k * w0**2) / (2*L)
sci = 0.5
    
pdf = theta(sci, sigma, 10)