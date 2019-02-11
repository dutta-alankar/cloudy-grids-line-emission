#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 18:26:56 2019

@author: alankar
"""

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('OIII_spectra.lin',skiprows=2)
#print(data)

incr = 0.05
ne = np.arange(np.min(data[:,0]),np.max(data[:,0])+incr,incr)
T = np.arange(np.min(data[:,1]),np.max(data[:,1])+incr,incr)

x, y, ratio = data[:,0], data[:,1], data[:,4]
counter = 0
nene, TT = np.meshgrid(ne,T)
ratio_cont = np.zeros((len(T),len(ne)))
for i in range(len(ne)):
    for j in range(len(T)):
        ratio_cont[j,i] = ratio[counter]
        counter += 1

plt.figure(figsize=(12,8))
# a base value of -30 is provided to avoid nans        
if (-30.0 in ratio): 
    import heapq
    plt.pcolor(ne,T,-ratio_cont,vmin=heapq.nsmallest(2,np.array(list(set(-ratio))))[-1])
    plt.colorbar()
    rcParams['contour.negative_linestyle'] = 'solid'
    CS = plt.contour(ne,T,-ratio_cont, 10,colors='white',vmin=heapq.nsmallest(2,np.array(list(set(-ratio))))[-1])
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
else: 
    plt.pcolor(ne,T,-ratio_cont) #change normalization with respect to 4363 A line intensity
    plt.colorbar()
    rcParams['contour.negative_linestyle'] = 'solid'
    CS = plt.contour(ne,T,-ratio_cont, 10,colors='white')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')

#plt.figure(figsize=(20,10))
plt.xlabel(r'$\log(n_H)$',size=18)
plt.ylabel(r'$\log(T)$',size=18)
plt.title(r'$\frac{I(\lambda = 5007 \AA)}{I(\lambda = 4363 \AA)}$',size=20,y=1.03)
plt.savefig('line_ratio.jpg')
plt.show()
