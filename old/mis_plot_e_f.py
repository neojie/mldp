#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:04:55 2020

@author: jiedeng
"""
import numpy as np
import matplotlib.pyplot as plt
out=np.loadtxt('std_vs_t.out')
plt.figure()
plt.plot(out[:,0],out[:,1],'o',label = 'e')
plt.plot(out[:,0],out[:,2],'^',label = 'f')
#plt.plot(out[:,0],out[:,3],'s',label = 'v')


plt.xlabel('sigma (eV)')
plt.ylabel('standard deviation')
plt.yscale('log')
plt.legend()
plt.show()

#np.savetxt('std_vs_t.out', out)
