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
X = np.load(filename_ind1 + '.npy')

filename_ind2 = "data/indeces_SCpspace_m_" + option 
Y = np.load(filename_ind2 + '.npy')


filename = "data/result_SCpspace_none"
baseline = np.load(filename + '.npy')



Y, X = np.meshgrid(Y, X)


fig = mlab.figure()



surf2 = mlab.surf(X, Y, data, colormap='Oranges', warp_scale="auto")
surf1 = mlab.surf(X, Y, baseline, colormap='Blues', warp_scale="auto")



#surf.actor.property.opacity = 0.5

#mlab.view(60, 74, 17, [-2.5, -4.6, -0.3])
#mlab.outline(surf1, color=(.7, .7, .7), extent=ax_extent)
mlab.axes(surf1, color=(1, 1, 1),xlabel='$\kappa$', ylabel='$\mu$', zlabel='z')


