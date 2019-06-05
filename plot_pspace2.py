# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter



def plot(option):
    
    filename = "data/result_pspace2_" + option 
    data = data2 = np.load(filename + '.npy')
    
    filename = "data/result_pspace2_p_" + option
    data2 = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_pspace2_k_" + option 
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = "data/indeces_pspace2_eta_" + option 
    etas = np.load(filename_ind2 + '.npy')
    
    filename = "data/result_pspace2_none"
    baseline = np.load(filename + '.npy')
    
    ks, etas = np.meshgrid(ks, etas)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(r"Transmitter side Scissors")
    ax.set_xlabel(r'$\kappa$', size=15)
    ax.set_ylabel(r'$\eta$', size=15)
    ax.set_zlabel(r'Key rate', size=10)
    #ax.set_zscale('log')
    
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
    
