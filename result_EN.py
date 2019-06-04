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
N = 20
option = 'rps'
eta = 0.01


## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

els = []
ps = []
ks = np.linspace(0000.1, .999, 20)
mus = np.linspace(0000.1, 1.5, 20)

for k in ks:
    el_temp = []
    p_temp = []
    for mu in mus:
        print("--->", k, mu)
        sys.reset_state(2)
        r = np.arcsinh(np.sqrt(mu))
        sys.apply_TMS(r, [0, 1])
                
    
        # Transmitter Scissors
        if option == 'tsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success)
        
        elif option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        elif option == 'none':
            p_success = 1
    
        elif option == 'tct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P SUCCESS:", p_success)
    
    
        sys.apply_loss_channel(eta, 1)
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success) 
            
        elif option == 'rps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
            
        elif option == 'rct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P SUCCESS:", p_success)

        el = measurements.log_neg(sys.state, [0, 1])
        print("Logarithmic Negativity:", el)
        el_temp += [el]
        p_temp += [p_success]
    els += [el_temp]
    ps += [p_temp]


# Save the resuls
filename = "data/result_ENpspace_" + option 
els = np.array(els)
np.save(filename, els)

filename = "data/result_ENpspace_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_ENpspace_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_ENpspace_m_" + option 
np.save(filename_ind2, mus)


############################################ PLOT


#option = 'tsc'

filename = "data/result_ENpspace_" + option 
data = data2 = np.load(filename + '.npy')

filename = "data/result_ENpspace_p_" + option
data2 = np.load(filename + '.npy')

filename_ind1 = "data/indeces_ENpspace_k_" + option 
ks = np.load(filename_ind1 + '.npy')

filename_ind2 = "data/indeces_ENpspace_m_" + option 
mus = np.load(filename_ind2 + '.npy')

#ks, mus = np.meshgrid(ks, mus)
mus, ks = np.meshgrid(mus, ks)


#fig = mlab.figure()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title(r"Logarithmic Negativity")
ax.set_xlabel(r'$\kappa$', size=15)
ax.set_ylabel(r'$\mu$', size=15)
ax.set_zlabel(r'$E_N$', size=10)


# Plot the surface.
#surf = mlab.surf(X, Y, data, colormap='Blues')
#surf.actor.property.opacity = 0.5

surf = ax.plot_surface(ks, mus, data, cmap=cm.coolwarm,
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