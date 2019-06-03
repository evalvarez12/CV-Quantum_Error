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
#from mayavi import mlab

############################################ CALCULATIONS

## Parameters
N = 15
mpne = 0.001
mu = 1.2 
f = 0.95
option = 'tsc'

r = np.arcsinh(np.sqrt(mu))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])
sys.save_state()

key_rates = []
ps = []
etas = np.logspace(-3, -1, base=10, num=20)
ks = np.linspace(0000.1, .999, 20)

for k in ks:
    k_temp = []
    p_temp = []
    for eta in etas:
        print("--->", k, mu)
        sys.load_state()
        
    
        # Transmitter Scissors
        if option == 'tsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success)
        
        if option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        if option == 'none':
            p_success = 1
    
        sys.add_TMSV(r_eve)

        theta = np.arccos(np.sqrt(eta))
        sys.apply_BS(theta, [1, 2])
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success) 
            
            
        if option == 'rps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)

        kr = measurements.key_rate(sys, f, p_success)
        print("Key rate:", kr)
        k_temp += [kr]
        p_temp += [p_success]
    key_rates += [k_temp]
    ps += [p_temp]


# Save the resuls
filename = "data/result_SCpspace2_" + option 
key_rates = np.array(key_rates)
np.save(filename, key_rates)

filename = "data/result_SCpspace2_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_SCpspace2_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_SCpspace2_eta_" + option 
np.save(filename_ind2, etas)


############################################ PLOT


#option = 'tsc'

filename = "data/result_SCpspace2_" + option 
data = data2 = np.load(filename + '.npy')

filename = "data/result_SCpspace2_p_" + option
data2 = np.load(filename + '.npy')

filename_ind1 = "data/indeces_SCpspace2_k_" + option 
ks = np.load(filename_ind1 + '.npy')

filename_ind2 = "data/indeces_SCpspace2_eta_" + option 
etas = np.load(filename_ind2 + '.npy')


filename = "data/result_SCpspace2_none"
baseline = np.load(filename + '.npy')



#ks, mus = np.meshgrid(ks, mus)
etas, ks = np.meshgrid(etas, ks)


#fig = mlab.figure()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title(r"Transmitter side Scissors")
ax.set_xlabel(r'$\kappa$', size=15)
ax.set_ylabel(r'$\eta$', size=15)
ax.set_zlabel(r'Key rate', size=10)
#ax.set_zscale('log')

# Plot the surface.
#surf = mlab.surf(X, Y, data, colormap='Blues')
#surf.actor.property.opacity = 0.5

#data = np.log(data)

surf = ax.plot_surface(ks, etas, data, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

#surf2 = ax.plot_surface(X, Y, data2,
#                       linewidth=0, antialiased=False)
#ax.set_zscale('log')


# Customize the z axis.a
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

#v1_matplotlib()
#v2_mayavi(False)
#v2_mayavi(True)