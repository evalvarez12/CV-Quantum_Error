# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:51:31 2019

@author: Eduardo Villasenor 
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

def plot(option):
    
    filename = "data/result_pspace_" + option 
    data = data2 = np.load(filename + '.npy')
    
    filename = "data/result_pspace_p_" + option
    data2 = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_pspace_k_" + option 
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = "data/indeces_pspace_m_" + option 
    mus = np.load(filename_ind2 + '.npy')
    
    ks, mus = np.meshgrid(ks, mus)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(r"Transmitter side Scissors")
    ax.set_xlabel(r'$\kappa$', size=15)
    ax.set_ylabel(r'$\mu$', size=15)
    ax.set_zlabel(r'Key rate', size=10)
    
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
    