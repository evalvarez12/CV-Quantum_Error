# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:15:44 2021

@author: z5239621
"""

import numpy as np

V = 50
T = 1
eps = 0 
g = 1
eta = 1
r = 0


A = (V + T*(V-1) + 1 + eps)*(g*eta)**2 + g**2*(1-eta**2)
B = -2 *g*eta*np.sqrt(T*(V**2-1))
C = (g*eta)**2 +1

F = 2/np.sqrt((A+B+C*(np.cosh(r) + np.sinh(r))**2)*(A+B+C*(np.cosh(r) - np.sinh(r))**2))
print(F)


print(2/(2*V - 2*np.sqrt(V**2-1) + 2))