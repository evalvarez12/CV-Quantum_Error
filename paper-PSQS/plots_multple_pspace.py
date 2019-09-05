# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:44:10 2019

@author: Eduardo Villaenor Alvarez
"""

import numpy as np
from mayavi import mlab
import src.names as names



def safe_log(x, minval=0.000000001):
    res_array = -np.log(x.clip(min=minval))
    res_array[x >= 1] = 0
    return res_array



def plotKR(option, N, eta):
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "KR"
    measurementp = "KR_p"

    
    # Save the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    datap = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    X = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    Y = np.load(filename_ind2 + '.npy')
    
    

    filename_base =  names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol='none')
    baseline = np.load(filename_base + '.npy')
   
    print("k:", min(X), max(X))
    print("mu:", min(Y), max(Y))
    xmin = min(X)
    xmax = max(X)
    xmin = 0 
    xmax = 1
    
    ymin = min(Y)
    ymax = max(Y)

    zmin = 0
    zmax = .000008
    
    factorz = 100000
    data = factorz*data
    
    ax_ranges = [xmin, xmax, ymin, ymax, zmin, zmax]
    ax_scale = [1, 1, factorz]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)
    
    Y, X = np.meshgrid(Y, X)
#    X, Y = np.mgrid[-1:1:20j, -1:1:20j]

    fig = mlab.figure()

    zero = np.zeros_like(data)

    surf1 = mlab.surf(X, Y, data, colormap='Oranges')
#    surf2 = mlab.surf(X, Y, 10*baseline, colormap='Blues')
    surf3 = mlab.surf(X, Y, zero, colormap='bone', opacity=.5)

    mlab.title("Photon Subtraction")
    mlab.outline(surf1, color=(.7, .7, .7), extent=ax_extent)
    mlab.axes(surf1, color=(.7, .7, .7), extent=ax_extent,
              ranges=ax_ranges,
              xlabel='$\kappa$', ylabel='$\mu$', zlabel='Key rate')

    surf1.actor.property.opacity = 1
    surf3.actor.property.opacity = 0.5
    fig.scene.renderer.use_depth_peeling = 1



#    ax_ranges = [-2, 2, -2, 2, 0, 8]
#    ax_scale = [1.0, 1.0, 0.4]
#    ax_extent = ax_ranges * np.repeat(ax_scale, 2)
#
#    surf3 = mlab.surf(mx, my, mz1, colormap='Blues')
#    surf4 = mlab.surf(mx, my, mz2, colormap='Oranges')
#
#    surf3.actor.actor.scale = ax_scale
#    surf4.actor.actor.scale = ax_scale
#    mlab.view(60, 74, 17, [-2.5, -4.6, -0.3])
#    mlab.outline(surf3, color=(.7, .7, .7), extent=ax_extent)
#    mlab.axes(surf3, color=(.7, .7, .7), extent=ax_extent,
#              ranges=ax_ranges,
#              xlabel='x', ylabel='y', zlabel='z')
#
#    if transparency:
#        surf3.actor.property.opacity = 0.5
#        surf4.actor.property.opacity = 0.5
#        fig.scene.renderer.use_depth_peeling = 1







def plotKR2(option, N, mu):
    k_name = 'var'
    mu_name = mu
    eta_name = 'var'
    measurement = "KR"
    measurementp = "KR_p"
    
    # Save the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    datap = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    X = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='eta')
    Y = np.load(filename_ind2 + '.npy')
    
    

    filename_base =  names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol='none')
    baseline = np.load(filename_base + '.npy')
   
    xmin = min(X)
    xmax = max(X)
    ymin = min(Y)
    ymax = max(Y)
    zmin = 0
    zmax = .08
    
    factorz = 10000000
    data = factorz*data
    
    ax_ranges = [xmin, xmax, ymin, ymax, zmin, zmax]
    ax_scale = [1, 1, factorz]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)
    
    
    zero = np.zeros_like(data)
    
    Y, X = np.meshgrid(Y, X)
#    X, Y = np.mgrid[-1:1:20j, -1:1:20j]

    fig = mlab.figure()
#
    
#    data = np.log(data)
#    baseline = np.log(baseline)
    
    surf2 = mlab.surf(X, Y, data, colormap='Oranges')
    surf1 = mlab.surf(X, Y, baseline, colormap='Blues')
    surf3 = mlab.surf(X, Y, zero, colormap='bone', opacity=1)
#    
    
    
    mlab.outline(surf1, color=(.7, .7, .7), extent=ax_extent)
    mlab.axes(surf1, color=(.7, .7, .7), extent=ax_extent,
              ranges=ax_ranges,
              xlabel='$\kappa$', ylabel='$\mu$', zlabel='Key rate')    
    
    surf1.actor.property.opacity = 1
    surf3.actor.property.opacity = 0.5
    fig.scene.renderer.use_depth_peeling = 1

    
def plotGRCI(option, N, eta):
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "GRCI"
    measurementp = "GRCI_p"
    
    # Save the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    datap = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    X = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    Y = np.load(filename_ind2 + '.npy')
    
    data[data < 0] = 0


    print("k:", min(X), max(X))
    print("mu:", min(Y), max(Y))
    xmin = min(X)
    xmax = max(X)
    ymin = min(Y)
    ymax = max(Y)
    zmin = 0
    zmax = .000008
    
    factorz = 100000
    data = factorz*data
    
    ax_ranges = [xmin, xmax, ymin, ymax, zmin, zmax]
    ax_scale = [1, 1, factorz]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)


#    filename_base =  names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol='none')
#    baseline = np.load(filename_base + '.npy')
   
    Y, X = np.meshgrid(Y, X)
#    X, Y = np.mgrid[-1:1:20j, -1:1:20j]


    fig = mlab.figure()

    surf2 = mlab.surf(X, Y, data, colormap='Oranges')
#    surf1 = mlab.surf(X, Y, baseline, colormap='Blues', warp_scale="auto")

    mlab.outline(surf1, color=(.7, .7, .7), extent=ax_extent)
    mlab.axes(surf1, color=(.7, .7, .7), extent=ax_extent,
              ranges=ax_ranges,
              xlabel='$\kappa$', ylabel='$\mu$', zlabel='Key rate')    
    
#    surf1.actor.property.opacity = 1
#    surf3.actor.property.opacity = 0.5
#    fig.scene.renderer.use_depth_peeling = 1
    
    
def plotEN(option, N, eta):
    k_name = 'var'
    mu_name = 'var'
    eta_name = eta
    measurement = "EN"
    measurementp = "EN_p"
    
    # Save the resuls
    filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
    data = np.load(filename + '.npy')
    
    filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
    datap = np.load(filenamep + '.npy')
    
    filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
    X = np.load(filename_ind1 + '.npy')
    
    filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
    Y = np.load(filename_ind2 + '.npy')
    
    

    filename_base =  names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol='none')
    baseline = np.load(filename_base + '.npy')
    
    filename_base2 =  names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol='rps')
    baseline2 = np.load(filename_base + '.npy')
   
   
    
    print("k:", min(X), max(X))
    print("mu:", min(Y), max(Y))
    xmin = min(X)
    xmax = max(X)
    xmin = 0 
    xmax = 1
    
    ymin = min(Y)
    ymax = max(Y)
    
    zmin = 0
    zmax = .17
    
    factorz = 2
    data = factorz * data
    baseline = factorz * baseline
    baseline2 = factorz * baseline2
    
    
    ax_ranges = [xmin, xmax, ymin, ymax, zmin, zmax]
    ax_scale = [1, 1, factorz]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)
    
    Y, X = np.meshgrid(Y, X)
#    X, Y = np.mgrid[-1:1:20j, -1:1:20j]
    
    fig = mlab.figure()

    surf1 = mlab.surf(X, Y, data, colormap='Oranges')
#    surf2 = mlab.surf(X, Y, baseline, colormap='Greys')
    surf3 = mlab.surf(X, Y, baseline2, colormap='Greens')


    mlab.outline(surf1, color=(.7, .7, .7), extent=ax_extent)
    mlab.axes(surf1, color=(.7, .7, .7), extent=ax_extent, ranges=ax_ranges,
              xlabel='$\kappa$', ylabel='$\mu$', zlabel='$E_N$')    
    
#    surf1.actor.property.opacity = 1
#    surf2.actor.property.opacity = 0.6
#    fig.scene.renderer.use_depth_peeling = 1
    
    
#plotKR('rps', N=20, eta=0.01)
#plotKR2('rps', N=20, mu=1.3)
#plotGRCI('rsc', N=20, eta=0.01)


plotEN('rsc', N=20, eta=0.01)