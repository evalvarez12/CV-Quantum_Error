# -*- coding: utf-8 -*-
"""
Function to plot 3d Wigner functions
Created on Fri Jul 12 14:22:24 2019

@author: Eduardo Villasenor
"""

# 3D Wigner functions for a series of cat states. #
import qutip as qt
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def plot(psi, borders, num=300):

    start, end = borders
    xvec = np.linspace(start, end, num)

    X, Y = np.meshgrid(xvec, xvec)

    W = qt.wigner(psi, xvec, xvec)

    fig = plt.figure()
    ax = Axes3D(fig, azim=-62.5, elev=25)

    ax.plot_surface(X, Y, W, rstride=2, cstride=2, cmap=cm.jet, lw=.1)
    #ax.set_xlim3d(-9,9)
    #ax.set_xlim3d(-9,9)
    #ax.set_zlim3d(-.2,0.2)
#    plt.margins(0,0)
    ax.set_axis_off()
    plt.show()
#    ax.set_frame_on('false')
    #    title(r'$| \psi >= \frac{1}{\sqrt{2}}(|\alpha>-|-\alpha>)$'+r' $\alpha=$'+str(round(x,2)))
    #    savefig("cat_state_"+str(k)+".png")


#setup constants:
N = 20  # size of the Hilbert space
a = qt.destroy(N)  # annihilation operator
alpha = 2.4
psi = (qt.displace(N, alpha) * qt.basis(N, 0)
       - qt.displace(N, -1*alpha) * qt.basis(N, 0)).unit()

# psi = qt.displace(N, alpha) * qt.basis(N, 0)
# psi = qt.squeeze(N, alpha) * qt.basis(N, 0)
# psi = qt.basis(N, 1)

plot(psi, [-6, 6])
