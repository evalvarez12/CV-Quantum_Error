# -*- coding: utf-8 -*-
"""
Code to test CV system
Created on Mon Apr 15 10:47:50 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import tools
import matplotlib.pyplot as plt
import theory



# Parameters
N = 30
mpn = 1
t = np.pi/3

wx = np.linspace(-5, 5)


# Testing two mode squeeze

r = np.arcsinh(np.sqrt(mpn))
print("Squeezing:", r)

## Initialize state
sys = cv.System(N, Nmodes=2, cm=True)
print("Initial CM:", sys.cm)


sys.apply_TMS(mpn, [0, 1])



print("Sys CM:", sys.cm)



sys.set_quadratures_basis()
CM = sys.get_full_CM()
print("CM:", CM)


S_theory = theory.two_mode_squeeze(r)
print("Theory S:", S_theory)

print("Theory CM:", tools.matrix_sandwich(S_theory, np.eye(4)))
