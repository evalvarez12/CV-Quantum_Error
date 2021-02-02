# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:18:53 2021

@author: z5239621
"""

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

def fidelity(V, T, tsc, eps, eta, g, alpha):
    gsc = np.sqrt((1 - tsc)/tsc)
    g = g*eta
    gt = g - 1
    
#    print(gsc, gt)
    
    A = 1/2 * (T*(V-1) + eps + 2)
    A2 = 1/ ((1 + gsc) *A)
    A3 = T * (V**2 - 1)/ (4*A)
    B1 = (gsc / (2 * A)) * np.sqrt(T * (V**2 - 1))
    A4 = 1/2 * (V + 2*A3 + g**2)
    A5 = 1 - gsc**2 * (1 - 1/A)
    A6 = -B1*g + gsc**2 * (A3 - g * (1 - 1/A))
    A7 = (gsc * g)**2 * A3
    A8 = g**2/2 - A4
    A9 = A8 + 1/2
    A10 = gt**2/A9
    
    
    F = 1/A2*np.exp(-A10 * np.abs(alpha)**2) * (A5/A9 + A6/(np.sqrt(2)*A9**3) \
        * (4*A9 - 4*gt**2*np.abs(alpha)**2) + A7/(2*A9**5) \
        * (28*A9**2 - 56*A9*(gt*np.abs(alpha))**2 + 16*gt**4*(np.real(alpha)**4 \
           + np.imag(alpha)**4 + (np.real(alpha)*np.imag(alpha))**2)))       
    return F
    
    
V = np.random.rand() + 1
T = np.random.rand()
tsc = np.random.rand() * 0.5
tsc = 0.2
T = 1

print(V, T, tsc)
print(fidelity(V, T, tsc, 0, 1, 1, 1+1j))