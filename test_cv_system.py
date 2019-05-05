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
import unittest

N = 10
sys = cv.System(N, Nmodes=1)
r = .6
sys.apply_SMD(r)

statei = sys.state
print(statei)
k = .1

g = np.sqrt(1/k - 1)
sys.apply_scissor_exact(k)
print(sys.state)

print(g)
#class TestCVSystemsMethods(unittest.TestCase):
#
#    def test_apply_scissor_exact(self):
#        kappa = .4
#        sys.apply_scissor_exact(kappa)
#        print(sys.state)
#
#
#if __name__ == '__main__':
#    unittest.main()
