# -*- coding: utf-8 -*-
"""
Object to simulate the atmospheric effects over CV states
Created on Thu Sep  5 14:37:49 2019

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.special as sp

class EllipticalModel:
    def __init__(self, lamb, w0, r0, scint, L):
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
        mean = np.log((1+2.96 * a)**2 / (self.sigma**2 * np.sqrt((1 + 2.96 * a)**2 + 1.2 * a)))
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
        
    
    def dist_phi(self, size):
        return np.random.uniform(0, np.pi, size)

    
    def get_dist_vals(self, size):
        thetas = self.dist_thetas(size)
        theta1 = thetas[:, 0]
        theta2 = thetas[:, 1]
        x = self.dist_pos(size)
        y = self.dist_pos(size)
        phi = self.dist_phi(size)
        return [theta1, theta2, x, y, phi]


    def W_eff(self, w1, w2, phi):
        r0 = self.r0
        
        W = sp.lambertw
        a = 4 * r0**2/ (w1 * w2)
        b = np.exp((r0/w1)**2 * (1+ 2 * np.cos(phi)**2)) 
        c = np.exp((r0/w2)**2 * (1+ 2 * np.sin(phi)**2))
        
        return np.sqrt(4*r0**2 / W(a * b * c))
    
    
    def T0(self, w1, w2):
        r0 = self.r0**2

        I0 = lambda z: np.real(sp.iv(0, z))
        a = I0(r0 * (1/(w1**2) - 1/(w2**2)))
        b = np.exp(-r0 * (1/(w1**2) + 1/(w2**2)))
        c = -2*(1 - np.exp(-(r0/2)*(1/w1 - 1/w2)**2))
        d = np.exp(-(((w1+w2)**2/np.abs(w1**2-w2**2))/self.R(1/w1 - 1/w2))**self.lamb(1/w1 - 1/w2))
        return 1 - a * b - c * d
    
    
    def R(self, w):
        I0 = lambda z: np.real(sp.iv(0, z))
        r0 = self.r0**2
        
        a = 1 - np.exp(-.5 * r0 * w**2)
        b = 1 - np.exp(-r0 * w**2) * I0(r0 * w**2)
        return np.log(2*a/b)**(-1/self.lamb(w))


    def lamb(self, w):
        I0 = lambda z: np.real(sp.iv(0, z))
        I1 = lambda z: np.real(sp.iv(1, z))
        r0 = self.r0**2

        a = 2 * r0 * w**2 * (np.exp(-r0 * w**2) * I1(r0 * w**2))              
        b = 1 - np.exp(-r0 * w**2) * I0(r0 * w**2)
        c = np.log((2*(1 - np.exp(-.5 * r0 * w**2)))/(1 - np.exp(-r0 * w**2) * I0(r0 * w**2)))
        return a * b / c


    def Te(self, w1, w2, x, y, phi):
        t0 = self.T0(w1, w2)
        phi0 = np.arctan(x/y)
        
        a = (np.sqrt(x**2 + y**2)/self.r0)/self.R(2/self.W_eff(w1, w2, phi - phi0))
        return t0 * np.exp(-a**self.lamb(2/self.W_eff(w1, w2, phi - phi0)))
    
    
    def get_te(self, size):
        theta1, theta2, x, y, phi = self.get_dist_vals(size)
        w1 = self.w0 * np.exp(theta1/2)
        w2 = self.w0 * np.exp(theta2/2)
        
        return np.real(self.Te(w1, w2, x, y, phi))
    
# All in mm
lamb = 1550 * 1e-6 # nm to mm
k = 2 * np.pi/(lamb)
w0 = 80
r0 = 110
L = 100 * 1e6 # km to mm
scint = 0.9

sigma = (k * w0**2) / (2*L)
print("Sigma:", sigma)

model = EllipticalModel(lamb, w0, a, scint, L)
    
ts = model.get_te(1000000)

import matplotlib.pyplot as plt

plt.hist(ts, 1000)
plt.show()


