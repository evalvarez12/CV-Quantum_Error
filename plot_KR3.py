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

options = ['none', 'tps', 'rps', 'tqs', 'rqs']

for option in options:
    
    
    if option == 'none':
        k =  0
    elif option[1:] == 'ps': 
        k = 0.95
    elif option == 'tqs':
        k = 0.05


    if option == 'rqs':
        N = 10

    params = ["r=" + str(0.93), "r_eve=" + str(0.033), "k=" + str(k)]

    filename = names.measurements_line(N, 'KR3', params, option)
    data = np.load(filename + ".npy")
    datas += [data]
    


none, tps, rps, tqs, rqs = datas

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
rqs[rqs < 0 ] = 0 
none[ none < 0] = 0 

#rps = rps - none
#tqs = tqs - none
#
#rps[rps < 0 ] = 0
#tqs[tqs < 0 ] = 0


rps[rps < none] = 0
tqs[tqs < none] = 0
rqs[rqs < none] = 0
rqs[rqs < tqs] = 0
none[none < tqs] = 0
none[none < rps] = 0

tps[tps < none] = 0


step_qs = 0.0005
step_ps = 0.0000005

none2[none2 > -.001] = 1
none2[none2 < -.001] = 0


none = np.log10(none)
rps = np.log10(rps)
tqs = np.log10(tqs)
r = np.log10(rqs)
#none2 = np.log10(none2)


#none2[rps > step_ps/1.2] = 0
#none2[tqs > step_qs/1.6] = 0
#none2[none2 != 0] = 1


rps = rps.reshape(50, 50)
tqs = tqs.reshape(50, 50)
none = none.reshape(50,50)
none2 = none2.reshape(50, 50)
tps = tps.reshape(50, 50)
rqs = rqs.reshape(50, 50)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_title( "tqs")
#surf = ax.plot_surface(y, x, none, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_title( "tps")
#surf = ax.plot_surface(y, x, tps, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title( "tqs")
surf = ax.plot_surface(y, x, tqs, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title( "rqs")
surf = ax.plot_surface(y, x, rqs, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)


fig = plt.figure()
pa = plt.imshow(rps, interpolation='nearest', cmap=cm.Reds)
cba = plt.colorbar(pa, shrink=1)
pb = plt.imshow(tqs, interpolation='nearest', cmap=cm.Greens)
cbb = plt.colorbar(pb, shrink=1)
cba.set_label('positive')
cbb.set_label('negative')

pc = plt.imshow(none, interpolation='nearest', cmap=cm.Greys)
cbc = plt.colorbar(pb, shrink=1)



fig = plt.figure()

m = np.amax(tqs)
levels1 = np.arange(0.0, m, step_qs) + step_qs


m = np.amax(rps)
levels2 = np.arange(0.0, m, step_ps) + step_ps


cmap_bw = colors.ListedColormap(['grey', 'lemonchiffon'])
bounds=[0,.5,1]
norm = colors.BoundaryNorm(bounds, cmap_bw.N)

cs = plt.contourf(x, y, none2, 1, cmap=cmap_bw, alpha=.5, norm=norm)

cs = plt.contour(x, y, none, 10, cmap='Greys', alpha=1)




cs = plt.contourf(x, y, rps, 30, cmap='Reds', alpha=1)
cs = plt.contourf(x, y, tqs, 30, cmap='Greens', alpha=1)
#cs = plt.contourf(x, y, rqs, 30, cmap='Greens', alpha=1)



plt.rcParams["font.family"] = "Times New Roman"


plt.text(5, 6.5, r'TMSV', size=13)
plt.text(0.6, 17, r'No key rate', size=13)
plt.text(.9, 3, r'tQS', size=13)
plt.text(6.3, 17, r'rPS', size=13)


cba = plt.colorbar(pa, shrink=.7)
cbb = plt.colorbar(pb, shrink=.7)
cbc = plt.colorbar(pc, shrink=.5)



tick_locator1 = ticker.MaxNLocator(nbins=5)
tick_locator2 = ticker.MaxNLocator(nbins=5)
tick_locator3 = ticker.MaxNLocator(nbins=6)

cba.locator = tick_locator1
cba.update_ticks()

cbb.locator = tick_locator2
cbb.update_ticks()

cbc.locator = tick_locator3
cbc.update_ticks()

cba.ax.tick_params(labelsize=13) 
cbb.ax.tick_params(labelsize=13) 
cbc.ax.tick_params(labelsize=13) 


plt.ylabel(r'Attenuation (dB)', size=13)
plt.xlabel(r'Squeezing (dB)', size=13)
plt.tick_params(axis='both', which='major', labelsize=12)


plt.show()






