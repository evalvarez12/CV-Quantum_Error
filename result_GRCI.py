# -*- coding: utf-8 -*-
"""
Replicating Gaussian RCI plots

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import measurements
import numpy as np
import scipy.io
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import copy

# Initial parameters
N = 10

# Sweep parameters
kappa = 0.005
mu = 0.05

# Fixed parameters
mu_aux = 0.01
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01


sys = cv.System(N, Nmodes=2, cm=False)
sys.save_state()
results = []
results2 = []

for kappa in np.linspace(0.0001, 0.01, 20):
    res = []
    res2 = []
    for mu in np.linspace(0.0001, .1, 20):
        
        r = np.arcsinh(np.sqrt(mu))
        sys.apply_TMS(r, pos=[0, 1])
        
        
        sys.apply_loss_channel(eta, 1)
        
        sys.apply_scissors(kappa, r_aux, 1)
        sys2 = copy.deepcopy(sys)
        
        rci = measurements.RCI(sys, [0])
        ci = measurements.RCI(sys2, [1])

        print(mu, kappa, rci, ci)
        
        sys.load_state()
        res += [rci]
        res2 += [ci]
    results += [res]
    results2 += [res2]


filename = "data/rci_plot_1NLA"

results = np.array(results)
np.save(filename, results)


filename2 = "data/rci_plot_1NLA2"

results2 = np.array(results2)
np.save(filename2, results2)

#################

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
#X = np.arange(-5, 5, 0.25)
X = np.linspace(0.0001, .1, 20)

#Y = np.arange(-5, 5, 0.25)
Y = np.linspace(0.0001, 0.01, 20)

X, Y = np.meshgrid(X, Y)

filename = "data/rci_plot_1NLA.npy"
Z = np.load(filename)

filename2 = "data/rci_plot_1NLA2.npy"
Z2 = np.load(filename2)

eta = 0.01
floorval = -np.log2(1 - eta)
floor = floorval * np.ones_like(Z)

Z[Z < floorval] = 0
Z2[Z2 < floorval] = 0

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

surf = ax.plot_surface(X, Y, Z2, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

norm = surf.norm


#floor = ax.plot_surface(X, YZ, floor, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False, norm=norm)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_ylabel(r"$\kappa$")
ax.set_xlabel(r"$\mu_{aux}$")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.cax.colorbar(im)
#ax.cax.toggle_label(True)


plt.show()










