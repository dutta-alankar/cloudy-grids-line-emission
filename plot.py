#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 18:26:56 2019

@author: alankar
"""

import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
import numpy as np
#import corner

data = np.loadtxt('OIII_spectra.lin',skiprows=2)
#print(data)

incrne = 0.01
incrT = 0.01
ne = np.arange(np.min(data[:,0]),np.max(data[:,0])+incrne,incrne)
T = np.arange(np.min(data[:,1]),np.max(data[:,1])+incrT,incrT)

x, y, ratio = data[:,0], data[:,1], data[:,5]
counter = 0
nene, TT = np.meshgrid(ne,T)
ratio_cont = np.zeros((len(T),len(ne)))
for i in range(len(ne)):
    for j in range(len(T)):
        ratio_cont[j,i] = ratio[counter]
        counter += 1

plt.figure(figsize=(13,10))
# a base value of -30 is provided to avoid nans        
if (-30.0 in ratio): 
    import heapq
    plt.pcolormesh(ne,T,-ratio_cont,vmin=heapq.nsmallest(2,np.array(list(set(-ratio))))[-1])
    plt.colorbar()
    rcParams['contour.negative_linestyle'] = 'solid'
    CS = plt.contour(ne,T,-ratio_cont, 10,colors='white',vmin=heapq.nsmallest(2,np.array(list(set(-ratio))))[-1])
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
else: 
    plt.pcolormesh(ne,T,-ratio_cont) #change normalization with respect to 4363 A line intensity
    plt.colorbar()
    rcParams['contour.negative_linestyle'] = 'solid'
    CS = plt.contour(ne,T,-ratio_cont, 10,colors='white')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

#plt.figure(figsize=(20,10))
plt.xlabel(r'$\log(n_H)$',size=18)
plt.ylabel(r'$\log(T)$',size=18)
plt.title(r'$\log\left(\frac{I(\lambda = 5007 \AA)}{I(\lambda = 4363 \AA)}\right)$',size=20,y=1.03)
rc('xtick', labelsize=18) 
rc('ytick', labelsize=18)
plt.savefig('line_ratio.png')
plt.show()

N = len(ne)
ratioT = np.matmul(-ratio_cont,np.ones((N,1)))/N
plt.figure(figsize=(13,10))
plt.plot(T,ratioT)
plt.grid()
plt.xlabel(r'$\log(T)$',size=18)
plt.ylabel(r'$\log\left(\frac{I(\lambda = 5007 \AA)}{I(\lambda = 4363 \AA)}\right)$',size=20)
plt.title(r'Line ratio averaged over Hydrogen density', size=21)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=12)
plt.savefig('line_ratio_Temperature.png')
plt.show()

"""
samples = np.column_stack((data[:,0], data[:,1], data[:,4]))
figure = corner.corner(samples, labels=[r"$\log(n_H)$", r"$\log(T)$"],
                       #quantiles=[0.16, 0.5, 0.84],
                       #show_titles=True, title_kwargs={"fontsize": 12},
                       )
plt.show()
"""
cool = np.log10(np.abs((10.**data[:,3] - 10.**data[:,4]))/10**(data[:,0]*data[:,2]))
counter = 0
nene, TT = np.meshgrid(ne,T)
cool_cont = np.zeros((len(T),len(ne)))
for i in range(len(ne)):
    for j in range(len(T)):
        cool_cont[j,i] = cool[counter]
        counter += 1

plt.figure(figsize=(13,10))
plt.pcolormesh(ne,T,cool_cont) #change normalization with respect to 4363 A line intensity
plt.colorbar()
rcParams['contour.negative_linestyle'] = 'solid'
CS = plt.contour(ne,T,cool_cont, 10,colors='white')
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
plt.xlabel(r'$\log(n_H)$',size=18)
plt.ylabel(r'$\log(T)$',size=18)
plt.title(r'Cooling Curve (log scale) $\frac{\Lambda}{n_e n_H} [erg.cm^3.s^{-1}]$',size=20,y=1.01)
rc('xtick', labelsize=18) 
rc('ytick', labelsize=18)
plt.savefig('cooling.png')
plt.show()
plt.show()