# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:04:09 2019

@author: Eduardo Villasenor
"""

import src.names as names
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors, ticker


N= 20
datas = []

options = ['none', 'rps', 'tqs']

for option in options:
    
    
    if option == 'none':
        k =  0
    elif option == 'rps': 
        k = 0.95
    elif option == 'tqs':
        k = 0.05


    params = ["r=" + str(0.93), "r_eve=" + str(0.033), "k=" + str(k)]

    filename = names.measurements_line(N, 'KR3', params, option)
    data = np.load(filename + ".npy")
    datas += [data]
    


none, rps, tqs = datas

vals1 = np.zeros_like(none)
vals2 = np.zeros_like(none)

vals1 = rps - none
vals2 = tqs - none

vals1[vals1 < 0] = 0
vals2[vals2 < 0] = 0

vals1 = vals1.reshape(50, 50)
vals2 = vals2.reshape(50, 50)




filename_ind1 = names.indeces_line(N, 'KR3', params, option, 'eta')
filename_ind2 = names.indeces_line(N, 'KR3', params, option, 'r')
    
etas = np.load(filename_ind1 + '.npy')
rs = np.load(filename_ind2 + '.npy')

x = -10*np.log10(np.exp(-2*rs))
y = -10*np.log10(etas)

x, y = np.meshgrid(x, y)


none2 = none.copy()

rps[rps < 0 ] = 0
tqs[tqs < 0 ] = 0
none[ none < 0] = 0 

rps = rps - none
tqs = tqs - none

rps[rps < 0 ] = 0
tqs[tqs < 0 ] = 0


step_qs = 0.0005
step_ps = 0.0000005

none2[none2 > 0] = 1
none2[none2 < 0] = 0
#none2[rps > step_ps/1.2] = 0
#none2[tqs > step_qs/1.6] = 0
#none2[none2 != 0] = 1


rps = rps.reshape(50, 50)
tqs = tqs.reshape(50, 50)
none = none.reshape(50,50)
none2 = none2.reshape(50, 50)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_title( "rps")
#surf = ax.plot_surface(y, x, rps, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_title( "tqs")
#surf = ax.plot_surface(y, x, none, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)






fig = plt.figure()
pa = plt.imshow(rps, interpolation='nearest', cmap=cm.Reds)
cba = plt.colorbar(pa, shrink=1)
pb = plt.imshow(tqs, interpolation='nearest', cmap=cm.Greens)
cbb = plt.colorbar(pb, shrink=1)
cba.set_label('positive')
cbb.set_label('negative')



fig = plt.figure()

m = np.amax(tqs)
levels1 = np.arange(0.0, m, step_qs) + step_qs


m = np.amax(rps)
levels2 = np.arange(0.0, m, step_ps) + step_ps


cmap_bw = colors.ListedColormap(['grey', 'white'])
bounds=[0,.5,1]
norm = colors.BoundaryNorm(bounds, cmap_bw.N)

cs = plt.contourf(x, y, none2, 1, cmap=cmap_bw, alpha=1, norm=norm)

cs = plt.contourf(x, y, rps, levels2, cmap='Reds', alpha=.7)
cs = plt.contourf(x, y, tqs, levels1, cmap='Greens', alpha=.7)



plt.rcParams["font.family"] = "Times New Roman"


plt.text(5, 6.5, r'TMSV', size=13)
plt.text(0.6, 17, r'No key rate', size=13)
plt.text(.9, 3, r'tQS', size=13)
plt.text(6.3, 17, r'rPS', size=13)


cba = plt.colorbar(pa, shrink=1)
cbb = plt.colorbar(pb, shrink=1)



tick_locator1 = ticker.MaxNLocator(nbins=5)
tick_locator2 = ticker.MaxNLocator(nbins=5)

cba.locator = tick_locator1
cba.update_ticks()

cbb.locator = tick_locator2
cbb.update_ticks()

cba.ax.tick_params(labelsize=13) 
cbb.ax.tick_params(labelsize=13) 


plt.xlabel(r'Attenuation (dB)', size=13)
plt.ylabel(r'Squeezing (dB)', size=13)
plt.tick_params(axis='both', which='major', labelsize=12)


plt.show()






