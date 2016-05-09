# -*- coding: utf-8 -*-
"""
Created on Tue May 03 22:14:25 2016

@author: Alexis
"""

import numpy as np
from scipy.special import binom as binom
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib
from time import localtime
from scipy.integrate import quad, dblquad
txtfile=open("resolution.txt")
tempvar=txtfile.readlines()

'''resolution'''
m= int(tempvar[0].split()[1])
n= int(tempvar[1].split()[1])
i= int(tempvar[2].split()[1])

rho=np.fromfile('rho.in',dtype=np.float64)
rho_dim=np.shape(rho)[0]
phi=np.fromfile('phi.in',dtype=np.float64)
phi_dim=np.shape(phi)[0]
#lm=np.fromfile('lm.in',dtype=np.float64).reshape(n,rho_dim)
#Rlm=np.fromfile('Rlm.in',dtype=np.float64).reshape(n,rho_dim)
Alm=np.fromfile('Alm.in',dtype=np.float64).reshape(n,phi_dim,rho_dim)
intensity=np.fromfile('intensity.in',dtype=np.float64).reshape(phi_dim,rho_dim)
population=np.fromfile('population.in',dtype=np.float64).reshape(phi_dim,rho_dim)

'''plots'''
save=False #set True if i want to save files automatically
R, P = np.meshgrid(rho, phi)
X, Y = R*np.cos(P), R*np.sin(P)


fig, (ax0) = plt.subplots(ncols=1, figsize=(9, 7))
cs0=ax0.pcolormesh(X, Y, intensity,cmap=plt.get_cmap('viridis'))
plt.colorbar(cs0)
#cs1=ax1.pcolormesh(X, Y, test,cmap=pyplot.get_cmap('viridis'),vmin=np.min(test), vmax=np.max(test))
fig = plt.figure(figsize=(11,9), dpi=100)
ax = plt.subplot(111, projection='3d')
ax.plot_surface(X, Y, intensity,  color="red", rstride=3, cstride=4, alpha=0.2)

fig, (ax0) = plt.subplots(ncols=1, figsize=(7, 6))
cs0=ax0.pcolormesh(X, Y, population,cmap=plt.get_cmap('viridis'))
plt.colorbar(cs0)
#cs1=ax1.pcolormesh(X, Y, test,cmap=pyplot.get_cmap('viridis'),vmin=np.min(test), vmax=np.max(test))
fig = plt.figure(figsize=(11,9), dpi=100)
ax = plt.subplot(111, projection='3d')
ax.plot_surface(X, Y, population,  color="red", rstride=3, cstride=4, alpha=0.2)
plt.show()
