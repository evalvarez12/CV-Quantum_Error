# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:51:31 2019

@author: Eduardo Villasenor 
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import src.names as names



def KR(option, N, eta):
    # File name parameters
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "KR"
    measurementp = "KR_p"

    # Load the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = data2 = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    data2 = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    mus = np.load(filename_ind2 + '.npy')
    
    mus, ks = np.meshgrid(mus, ks)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(option + str(N))
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
    
    
def KR2(option, N, mu):
    # File name parameters
    k_name = 'var'
    mu_name = mu
    eta_name = 'var'
    measurement = "KR"
    measurementp = "KR_p"

    # Load the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = data2 = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    data2 = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='eta')
    etas = np.load(filename_ind2 + '.npy')
    
    
    
    etas, ks = np.meshgrid(etas, ks)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(option + str(N))
    ax.set_xlabel(r'$\kappa$', size=15)
    ax.set_ylabel(r'$\eta$', size=15)
    ax.set_zlabel(r'Key rate', size=10)
    
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
    
    
def EN(option, N, eta):
        # File name parameters
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "EN"
    measurementp = "EN_p"
    
    # Load the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = data2 = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    data2 = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    mus = np.load(filename_ind2 + '.npy')
    
    mus, ks = np.meshgrid(mus, ks)
      
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(option + str(N))
    ax.set_xlabel(r'$\kappa$', size=15)
    ax.set_ylabel(r'$\mu$', size=15)
    ax.set_zlabel(r'$E_N$', size=10)
    
    
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
    

def GRCI(option, N, eta):
    # File name parameters
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "GRCI"
    measurementp = "GRCI_p"
    
    # Load the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = data2 = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    data2 = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    ks = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    mus = np.load(filename_ind2 + '.npy')
    
    mus, ks = np.meshgrid(mus, ks)
      
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(option + str(N))
    ax.set_xlabel(r'$\kappa$', size=15)
    ax.set_ylabel(r'$\mu$', size=15)
    ax.set_zlabel(r'$GRCI$', size=10)
    
    
    data[data < 0] = 0
    
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
    





    
#KR('rsc', 10, 0.01)
#KR('rps', 10, 0.01)

KR2('rsc', 10, 0.01)
KR2('rsc', 10, 1.3)

#GRCI('none', 20, 0.01)
#EN('rps', 20, 0.01)
#EN('none', 20, 0.01)
#EN('rsc', 20, 0.01)