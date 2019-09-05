# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import cv_system as cv
import measurements
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


############################################ CALCULATIONS

## Parameters
N = 10
mpn = 1.3
mpne = 0.001
f = 0.95
option = 'tsc'
eta = 0.5


theta = np.arccos(np.sqrt(eta))
r = np.arcsinh(np.sqrt(mpn))
r_eve = np.arcsinh(np.sqrt(mpne))
r_aux = np.arcsinh(np.sqrt(m_aux))

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])

# Save current state of the system
sys.save_state()

key_rates = []
ps = []
ks = np.linspace(0.001, 1, 20)
m_auxs = np.linspace(0.001, 1, 20)

for k in ks:
    k_temp = []
    p_temp = []
    for m_aux in m_auxs:
        r_aux = np.arcsinh(np.sqrt(m_aux))
        print("--->", k, m_aux)
        
        sys.load_state()
    
        # Transmitter Scissors
        if option == 'tsc':
            p_success = sys.apply_scissors_inverted(k, r_aux, 1)
            print("P SUCCESS:", p_success)
        #    print(sys.state)
    
        sys.add_TMSV(r_eve)

    
        sys.apply_BS(theta, [1, 2])
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors(k, r_aux, 1)
            print("P SUCCESS:", p_success) 
            

        kr = measurements.key_rate_nosimple(sys, f, p_success)
        print("Key rate:", kr)
        k_temp += [kr]
        p_temp += [p_success]
    key_rates += [k_temp]
    ps += [p_temp]


# Save the resuls
filename = "data/result_SCpspace_" + option 
key_rates = np.array(key_rates)
np.save(filename, key_rates)

filename = "data/result_SCpspace_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_SCpspace_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_SCpspace_m_" + option 
np.save(filename_ind2, m_auxs)


############################################ PLOT


option = 'tsc'

filename = "data/result_SCpspace_" + option 
data = data2 = np.load(filename + '.npy')

filename = "data/result_SCpspace_p_" + option
data2 = np.load(filename + '.npy')

filename_ind1 = "data/indeces_SCpspace_k_" + option 
X = np.load(filename_ind1 + '.npy')

filename_ind2 = "data/indeces_SCpspace_m_" + option 
Y = np.load(filename_ind2 + '.npy')

X, Y = np.meshgrid(X, Y)

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.set_title(r"Transmitter side Scissors")
ax.set_xlabel(r'$\kappa$', size=15)
ax.set_ylabel(r'$\mu_{aux}$', size=15)
ax.set_zlabel(r'Key rate', size=10)


# Plot the surface.
surf = ax.plot_surface(X, Y, data, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#surf2 = ax.plot_surface(X, Y, data2,
#                       linewidth=0, antialiased=False)
#ax.set_zscale('log')


# Customize the z axis.a
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
