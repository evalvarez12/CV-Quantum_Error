# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:59:19 2019

@author: Eduardo Villasenor
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


def plot():
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
    
    Z[Z < 0] = 0
    Z2[Z2 < 0] = 0
    
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
