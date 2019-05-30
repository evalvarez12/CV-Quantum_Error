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

################################## CALCULATIONS

# Initial parameters
N = 8

# Fixed parameters
mu_aux = 0.01
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01


sys1 = cv.System(N, Nmodes=2, cm=False)

sys2 = cv.System(N, Nmodes=2, cm=False)

sys1.save_state()
results_sc = []
results_ps = []
probabilities_sc = []
probabilities_ps = []

for kappa in np.linspace(0.0001, 0.01, 20):
#for kappa in np.linspace(0.0001, 0.03, 20):
    res_sc = []
    res_ps = []
    ps_sc = []
    ps_ps = []
    for mu in np.linspace(0.0001, .1, 20):
#    for mu in np.linspace(0.0001, .6, 20):
        
        r = np.arcsinh(np.sqrt(mu))
        sys1.apply_TMS(r, pos=[0, 1])
        
        
        sys1.apply_loss_channel(eta, 1)
        
#        sys2.set_state(sys1.state)
        
        
        # Quantum scissors 
#        p_sc = sys1.apply_scissors(kappa, r_aux, 1)
        p_sc = sys1.apply_scissors_exact(kappa, 1)

#        p_sc = sys1.apply_scissors_options(kappa, r_aux, 1, 'b')

        # Photon substraction
#        p_ps = sys2.apply_photon_subtraction(1- kappa, 1)
        p_ps = 1
        
#        print("Reversed")
        rci_sc = measurements.CI(sys1, [0])
#        print("Coherent")
#        ci = measurements.CI(sys2, [1])

#        rci_ps = measurements.CI(sys2, [0])
        rci_ps = 0

        print(mu, kappa, p_sc, rci_sc)
#        
#        print(rci > ci)
#        if rci > ci:
#            marker += 1
        
        sys1.load_state()
        
        res_sc += [rci_sc]
        res_ps += [rci_ps]
        ps_sc += [p_sc]
        ps_ps += [p_ps]
    
    results_sc += [res_sc]
    results_ps += [res_ps]
    probabilities_sc += [ps_sc]
    probabilities_ps += [ps_ps]

filename_sc = "data/rci_plot_SC"

results_sc = np.array(results_sc)
np.save(filename_sc, results_sc)

filename_p_sc = "data/p_plot_SC"

probabilities_sc = np.array(probabilities_sc)
np.save(filename_p_sc, probabilities_sc)



filename_ps = "data/rci_plot_PS"

results_ps = np.array(results_ps)
np.save(filename_ps, results_ps)

filename_p_ps = "data/p_plot_PS"

probabilities_ps = np.array(probabilities_ps)
np.save(filename_p_ps, probabilities_ps)

#################

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.linspace(0.0001, .1, 20)
#X = np.linspace(0.0001, .6, 20)

Y = np.linspace(0.0001, 0.01, 20)
#Y = np.linspace(0.0001, 0.03, 20)


X, Y = np.meshgrid(X, Y)

filename_sc = "data/rci_plot_SC.npy"
Z = np.load(filename_sc)

filename_ps = "data/rci_plot_PS.npy"
Z2 = np.load(filename_ps)

eta = 0.01
floorval = -np.log2(1 - eta)
floor = floorval * np.ones_like(Z)

Z[Z < floorval] = 0
Z2[Z2 < floorval] = 0

#Z = Z/2

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#surf = ax.plot_surface(X, Y, Z2, cmap=cm.bwr,
#                       linewidth=0, antialiased=False)

norm = surf.norm


#floor = ax.plot_surface(X, Y, floor, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False, norm=norm)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_ylabel(r"$\kappa$")
ax.set_xlabel(r"$\mu$")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.cax.colorbar(im)
#ax.cax.toggle_label(True)

############## FIGURE 2 with p_success

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
#X = np.arange(-5, 5, 0.25)
X = np.linspace(0.0001, .1, 20)

#Y = np.arange(-5, 5, 0.25)
Y = np.linspace(0.0001, 0.01, 20)

X, Y = np.meshgrid(X, Y)

filename_sc = "data/p_plot_SC.npy"
Z = np.load(filename_sc)

filename_ps = "data/p_plot_PS.npy"
Z2 = np.load(filename_ps)


# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#surf = ax.plot_surface(X, Y, Z2, cmap=cm.bwr,
#                       linewidth=0, antialiased=False)


# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_ylabel(r"$\kappa$")
ax.set_xlabel(r"$\mu$")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.cax.colorbar(im)
#ax.cax.toggle_label(True)



plt.show()










