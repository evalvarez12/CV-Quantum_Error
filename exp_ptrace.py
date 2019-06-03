# -*- coding: utf-8 -*-
"""
Example to plot nice Wigner functions

Created on Fri May 17 11:05:01 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import matplotlib.pyplot as plt
import measurements

N = 5
modes = 5

state = qt.tensor([qt.basis(N)]*modes)


traced = state.ptrace([1, 2])
