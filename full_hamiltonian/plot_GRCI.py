# -*- coding: utf-8 -*-
"""
Ploting routines to recreate RCI plots of quantum scissors paper

Created on Thu Mar  7 17:00:02 2019

@author: Eduardo Villasenor 
"""


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
#X = np.arange(-5, 5, 0.25)
X = np.linspace(0.001, 0.01, 20)

#Y = np.arange(-5, 5, 0.25)
Y = np.linspace(0.001, .1, 20)

X, Y = np.meshgrid(X, Y)

filename = "rci_plot_1NLA.npy"
Z = np.load(filename)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
