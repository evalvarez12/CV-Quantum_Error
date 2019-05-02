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


N = 10
a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))

class TestCVSystemsMethods(unittest.TestCase):
    

if __name__ == '__main__':
    unittest.main()