# -*- coding: utf-8 -*-
"""
Object to simulate the atmospheric effects over CV states
Created on Thu Sep  5 14:37:49 2019

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.special as sp

class EllipticalModel:
    def __init__(self, lamb, w0, r0, scint, L)
        self.k = 2 * np.pi / lamb
        self.w0 = w0
        self.r0 = r0
        self.scint = scint
        self.L = L
        self.sigma = (self.k * self.w0**2) / (2 * self.L) 
    
    
    def set_distance(self, L):
        self.L = L
        self.sigma = (self.k * self.w0**2) / (2 * self.L) 
    
    
    def dist_thetas(self, size):
        a  = self.scint**2 * self.sigma**(5/6)
        mean = np.log((1+2.96 * a)**2 / (sigma**2 * np.sqrt((1 + 2.96 * a)**2 + 1.2 * a)))
        var = np.log(1 + (1.2 * a) / (1 + 2.96 * a)**2)
        var12 = np.log(1 - (0.8 * a) / (1 + 2.96 * a)**2)
        
        mean_vec = [mean, mean]
        var_vec = np.array([[var, var12], [var12, var]])
        print(mean)
        print( var_vec)
        return np.random.multivariate_normal(mean_vec, var_vec, size)
        
        
    def dist_pos(self, size):
        a  = self.w0**2 * self.scint**2 * self.sigma**(-7/6)
        mean = 0
        var = 0.33 * a
        return np.random.normal(mean, var, size)
        
    
    def dist_phi(size):
        return np.random.uniform(0, np.pi, size)


    def W_eff(self, theta1, theta2, phi):
        w1 = self.w0 * np.exp(theta1/2)
        w2 = self.w0 * np.exp(theta2/2)
        
        W = sp.lambertw
        a = 4 * r0**2/ (w1 * w2)
        b = np.exp((r0/w1)**2 * (1+ 2 * np.cos(phi)**2)) 
        c = np.exp((r0/w2)**2 * (1+ 2 * np.sin(phi)**2))
        
        return 4*r0**2 / W(a * b * c)
    
    def T0()



k = 2 * np.pi/(1550 * 1e-6)
w0 = 80
a = 110
L = 1 * 1e6

sigma = (k * w0**2) / (2*L)
sci = 0.9
print("Sigma:", sigma)
pdf = thetas(sci, sigma, 10)