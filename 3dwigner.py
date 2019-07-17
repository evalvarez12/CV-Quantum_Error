# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 14:22:24 2019

@author: z5239621
"""

# 3D Wigner functions for a series of cat states. #   
from qutip.states import *   
from qutip.qobj import *   
from qutip.operators import *   
from qutip.wigner import *   
from mpl_toolkits.mplot3d import Axes3D   
from matplotlib import cm   
from pylab import *  
import qutip as qt

#setup constants: 
N = 20 # size of the Hilbert space 
#alpha=linspace(0.01,3,60) # values for the coherent state amplitude 
alpha = [2]

a = destroy(N) #annihilation operator 

#D = displace(N,1) # Displacement 


xvec = linspace(-9,9,300) 

X,Y = meshgrid(xvec, xvec)

k=1 # counter 


psi=basis(N,1).unit();
psi = qt.squeeze(N, .8)* basis(N)
#psi = qt.displace(N, 2)* basis(N)
# plot the wigner function 
W=wigner(psi,xvec,xvec) 
fig =figure() 
ax = Axes3D(fig,azim=-62.5,elev=25) 
ax.plot_surface(X, Y, W, rstride=2, cstride=2, cmap=cm.jet,lw=.1) 
ax.set_xlim3d(-9,9) 
ax.set_xlim3d(-9,9) 
ax.set_zlim3d(-.2,0.2) 
ax.set_axis_off() 
ax.set_frame_on('false') 
#    title(r'$| \psi >= \frac{1}{\sqrt{2}}(|\alpha>-|-\alpha>)$'+r' $\alpha=$'+str(round(x,2))) 
#    savefig("cat_state_"+str(k)+".png") 
k=k+1;