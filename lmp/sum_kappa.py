#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:39:51 2021

@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
ds = ['/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/4k/rp5/160-cpu/kappa/160/nvtnve',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/4k/r1/kappa/160/nvtnve_2',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/6k/rp5/kappa/160/nvtnve',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/4k/rp5/160-cpu/kappa/540/nvtnve',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/4k/rp5/160-cpu/kappa/1280/om8']
ps = [140, 3, 150]
ts = [4000, 4000, 6000]

ps = [140, 3, 150]
ts = ['4000 K', 4000, 6000]

fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)

kappa = np.loadtxt(ds[0] + '/log.kappa')
ax[0].plot(kappa[:,0], kappa[:,1],label= '4000 K, V/Vx=0.5')
ax[1].plot(kappa[:,0], kappa[:,2],label= '4000 K, V/Vx=0.5')
kappa = np.loadtxt(ds[1] + '/log.kappa')
ax[0].plot(kappa[:,0], kappa[:,1],label= '4000 K, V/Vx=1')
ax[1].plot(kappa[:,0], kappa[:,2],label= '4000 K, V/Vx=1')
kappa = np.loadtxt(ds[2] + '/log.kappa')
ax[0].plot(kappa[:,0], kappa[:,1],label= '6000 K, V/Vx=0.5')
ax[1].plot(kappa[:,0], kappa[:,2],label= '6000 K, V/Vx=0.5')
kappa = np.loadtxt(ds[3] + '/log.kappa')
ax[0].plot(kappa[:,0], kappa[:,1],label= '540, 4000 K, V/Vx=0.5')
ax[1].plot(kappa[:,0], kappa[:,2],label= '540, 4000 K, V/Vx=0.5')
kappa = np.loadtxt(ds[4] + '/log.kappa')
ax[0].plot(kappa[:,0], kappa[:,1],label= '1280, 4000 K, V/Vx=0.5')
ax[1].plot(kappa[:,0], kappa[:,2],label= '540, 4000 K, V/Vx=0.5')

ax[0].grid(True)
ax[0].legend()
ax[0].set_xscale('log')

ax[1].set_xlabel('dt (ps)')
ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
ax[1].grid(True)
ax[1].legend()
ax[1].set_xscale('log')
plt.show()

#for i in range(3):
#    kappa = np.loadtxt(ds[i] + '/log.kappa')
#    ax[0].plot(kappa[:,0], kappa[:,1],label= ps)
#    ax[1].plot(kappa[:,0], kappa[:,2],label='average')
#    ax[0].set_xlabel('dt (ps)')
#    ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
#ax[0].grid(True)
#ax[0].legend()
#ax[0].set_xscale('log')
#
#
#
#
#ax[1].plot(dt, cumsum_JxJx,label='x') # 
#ax[1].plot(dt, cumsum_JyJy,label='y')
#ax[1].plot(dt, cumsum_JzJz,label='z')
#ax[1].plot(dt, cumsum_JJ,label='average')
#
#ax[1].set_xlabel('dt (ps)')
#ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
#ax[1].grid(True)
#ax[1].legend()
#ax[1].set_xscale('log')
#plt.show()