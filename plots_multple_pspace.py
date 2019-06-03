# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:44:10 2019

@author: Eduardo Villaenor Alvarez
"""

import numpy as np
from mayavi import mlab


option = 'rps'

filename = "data/result_SCpspace_" + option 
data = data2 = np.load(filename + '.npy')

filename = "data/result_SCpspace_p_" + option
data2 = np.load(filename + '.npy')

filename_ind1 = "data/indeces_SCpspace_k_" + option 
ks = np.load(filename_ind1 + '.npy')

filename_ind2 = "data/indeces_SCpspace_m_" + option 
etas = np.load(filename_ind2 + '.npy')


filename = "data/result_SCpspace_none"
baseline = np.load(filename + '.npy')



etas, ks = np.meshgrid(etas, ks)


fig = mlab.figure()

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_title(r"Transmitter side Scissors")
#ax.set_xlabel(r'$\kappa$', size=15)
#ax.set_ylabel(r'$\eta$', size=15)
#ax.set_zlabel(r'Key rate', size=10)
#ax.set_zscale('log')


x = np.arange(-2, 2, 0.1)
y = np.arange(-2, 2, 0.1)
mx, my = np.meshgrid(x, y, indexing='ij')
mz1 = np.abs(mx) + np.abs(my)
mz2 = mx ** 2 + my ** 2


#data = np.ones_like(data)
#etas = np.round(etas, 2)
#ks = np.round(ks, 2)
#data = np.round(data, 2)
# Plot the surface.

surf2 = mlab.surf(ks, etas, data, colormap='Oranges', warp_scale="auto")
surf1 = mlab.surf(ks, etas, baseline, colormap='Blues', warp_scale="auto")


#surf = mlab.surf(mx, my, mz1, colormap='Blues', warp_scale="auto")

#surf.actor.property.opacity = 0.5

#ax_ranges = [-2, 2, -2, 2, 0, 8]
#ax_scale = [1.0, 1.0, 0.4]
#ax_extent = ax_ranges * np.repeat(ax_scale, 2)
#
#surf3 = mlab.surf(mx, my, mz1, colormap='Blues')
#surf4 = mlab.surf(mx, my, mz2, colormap='Oranges')
#
#surf3.actor.actor.scale = ax_scale
#surf4.actor.actor.scale = ax_scale
#mlab.view(60, 74, 17, [-2.5, -4.6, -0.3])
#mlab.outline(surf3, color=(.7, .7, .7), extent=ax_extent)
#mlab.axes(surf3, color=(.7, .7, .7), extent=ax_extent,
#          ranges=ax_ranges,
#          xlabel='x', ylabel='y', zlabel='z')





#data = np.log(data)

#surf = ax.plot_surface(ks, etas, data, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

#surf2 = ax.plot_surface(X, Y, data2,
#                       linewidth=0, antialiased=False)
#ax.set_zscale('log')


# Customize the z axis.a
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
#plt.show()
#mlab.show()