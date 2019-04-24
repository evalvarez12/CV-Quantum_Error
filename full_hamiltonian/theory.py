# -*- coding: utf-8 -*-
"""
Functions for symplectic representations and covariance matrices calculated by
hand for comparisson

Created on Wed Apr 24 16:39:28 2019

@author: Eduardo Villasenor
"""

import numpy as np


def phase_shift(theta):
    cc = np.cos(theta)
    ss = np.sin(theta)
    S =np.array([[cc, -ss, 0, 0], [ss, cc, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return S    


def beam_splitter(z):
    c1 = np.cos(np.absolute(z))
    s1 = np.sin(np.absolute(z))
    c2 = np.cos(np.angle(z))
    s2 = np.sin(np.angle(z))

    S =np.array([[c1, 0, s1*c2, s1*s2], [0, c1, -s1*s2, s1*c2], [-s1*c2, -s1*s2, c1, 0], [s1*s2, -s1*c2, 0, c1]])
    return S


def single_mode_squeeze(z):
    cc = np.cos(np.angle(z))
    ss = np.sin(np.angle(z))
    ch = np.cosh(np.absolute(z))
    sh = np.sinh(np.absolute(z))
    S =np.array([[ch + cc*sh, ss*sh, 0, 0], [ss*sh, ch- cc*sh, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return S


def two_mode_squeeze(z):
    cc = np.cos(np.angle(z))
    ss = np.sin(np.angle(z))
    ch = np.cosh(np.absolute(z))
    sh = np.sinh(np.absolute(z))
    S =np.array([[ch, 0, sh*cc, sh*ss], [0, ch, sh*ss, -sh*cc], [sh*cc, sh*ss, ch, 0], [sh*ss, -sh*cc, 0, ch]])
    return S