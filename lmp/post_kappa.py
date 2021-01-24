#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:39:51 2021

@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt

ds = ['/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/3k/r1/kappa/nvtnve_om8',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/4k/r1/kappa/160/nvtnve_om8',
      '/u/project/ESS/lstixrud/jd848/pv+hf/dp-train/lmp_run/6k/r1/kappa/nvtnve_om8']

VVx = [1, 1, 1]
T   = [3000, 4000, 6000]

fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
for i in range(len(ds)):
    kappa = np.loadtxt(ds[i] + '/log.kappa')
    label = str(VVx[i]) + ' ' + str(T[i]) + ' K' 
    ax[0].plot(kappa[:,0], kappa[:,1],label= label)
    ax[1].plot(kappa[:,0], kappa[:,2],label= label)
ax[0].grid(True)
ax[0].legend()
ax[0].set_xscale('log')

ax[1].set_xlabel('dt (ps)')
ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
ax[1].grid(True)
ax[1].legend()
ax[1].set_xscale('log')
plt.show()    

