# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:44:10 2019

@author: Eduardo Villaenor Alvarez
"""

import numpy as np
from mayavi import mlab
import src.names as names


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
   
#    Y, X = np.meshgrid(Y, X)
    X, Y = np.mgrid[-1:1:20j, -1:1:20j]

    fig = mlab.figure()

    zero = np.zeros_like(data)

    data = 100* data
    baseline = 10* baseline

    surf2 = mlab.surf(X, Y, data, colormap='Oranges')
    surf1 = mlab.surf(X, Y, baseline, colormap='Blues')
    surf3 = mlab.surf(X, Y, zero, colormap='bone', opacity=1)


    mlab.axes(surf1, color=(1, 1, 1),xlabel='$\kappa$', ylabel='$\mu$', zlabel='z')



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
   
    
    zero = np.zeros_like(data)
    
#    Y, X = np.meshgrid(Y, X)
    X, Y = np.mgrid[-1:1:20j, -1:1:20j]

    fig = mlab.figure()
#
    
#    data = np.log(data)
#    baseline = np.log(baseline)
    
    data = 100* data
    baseline = 10* baseline
    
    surf2 = mlab.surf(X, Y, data, colormap='Oranges')
    surf1 = mlab.surf(X, Y, baseline, colormap='Blues')
    surf3 = mlab.surf(X, Y, zero, colormap='bone', opacity=1)
#    
#    mlab.axes(surf1, color=(1, 1, 1),xlabel='$\kappa$', ylabel='$\eta$', zlabel='z')
    
    
    
def plotGRCI(option):
    k_name = 'var'
    mu_name = 'var'
    eta_name = 0.01
    measurement = "GRCI"
    measurementp = "GRCI_p"
    N = 20
    
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
   
    Y, X = np.meshgrid(Y, X)

    fig = mlab.figure()

    surf2 = mlab.surf(X, Y, data, colormap='Oranges', warp_scale="auto")
    surf1 = mlab.surf(X, Y, baseline, colormap='Blues', warp_scale="auto")

    mlab.axes(surf1, color=(1, 1, 1),xlabel='$\kappa$', ylabel='$\mu$', zlabel='z')
    
    
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
   
#    Y, X = np.meshgrid(Y, X)
    X, Y = np.mgrid[-1:1:20j, -1:1:20j]
    
    fig = mlab.figure()

    surf2 = mlab.surf(X, Y, data, colormap='Oranges')
    surf1 = mlab.surf(X, Y, baseline, colormap='Blues')

    mlab.axes(surf1, color=(1, 1, 1),xlabel='$\kappa$', ylabel='$\mu$', zlabel='z')
    
    
    
    
#plotKR('rps', N=20, eta=0.01)
#plotKR2('rps', N=20, mu=1.3)
plotEN('rsc', N=20, eta=0.01)