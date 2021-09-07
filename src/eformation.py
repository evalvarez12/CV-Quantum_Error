# -*- coding: utf-8 -*-
"""
Fucntions to compute the entanglement of formation of a Gaussian CV state

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.linalg as la
import qutip as qt


def GEOF(cm):
    """
    Compute the Gaussian entanglement of formation
    Covariance matrix =
    (a  0   x   0)
    (0  b   0   y)
    (x  0   c   0)
    (0  y   0   d)

    Ref: Quantum Inf Process (2015) 14:4179–4199
    """
    E = f(Delta(cm))
    return E

def Delta(cm):
    a, b, c, d, x, y = cm
    """
    Returns Delta prime
    """
    if ((c*b - a*d) == 0) and ((-a-d+c+b) == 0):
        x0 = 1
    else:
        x0 = (c*b - a*d) / (-a -d + c + b)

    a02 = (c-x0)/(a-x0)
    a022 = (d-x0)/(b-x0)
    print(a02, a022)

    b0 = np.sqrt(1 - 4/(a02 + 1/a02)**2)

    D = (a02*(a+b)/2 + (c+d)/(2*a02) - (x - y))/(a02 + 1/a02)
    Dp = (D + np.sqrt(D**2 - b0**2)) / (1 + np.sqrt(1-b0**2))
    return Dp

def f(x):
    cp = (1/np.sqrt(x) + np.sqrt(x))**2/4
    cm = (1/np.sqrt(x) - np.sqrt(x))**2/4

    return cp * np.log(cp) - cm * np.log(cm)


r1 = 1.2
r2 = 1.2

test_values = np.array([[2, 1.5, 1.2, -1],
                        [2, 1.5, 1, -1],
                        [3, 2, 1.8, -1.2],
                        [2.6, 1.7, 1.3, -0.9],
                        [3,2, 1.7, -1.2],
                        [2.5, 2, 1.3, -1.2]])

test_cm = np.zeros((6, 6))
test_cm[:, 0] = test_values[:, 0]*r1
test_cm[:, 1] = test_values[:, 0]/r1
test_cm[:, 2] = test_values[:, 1]*r2
test_cm[:, 3] = test_values[:, 1]/r2
test_cm[:, 4] = test_values[:, 2]*np.sqrt(r1*r2)
test_cm[:, 5] = test_values[:, 3]/np.sqrt(r1*r2)

print(test_cm)

print(GEOF(test_cm[0]))
print(GEOF(test_cm[1]))
print(GEOF(test_cm[2]))
print(GEOF(test_cm[3]))
print(GEOF(test_cm[4]))
print(GEOF(test_cm[5]))
