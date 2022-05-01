#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 21:20:52 2022

@author: jiedeng
"""
import sys
try:
    sys.path.insert(1, '/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/mldp/similiarity/')
except:
    print("run in local")

from stat_lib import show,select_chi,select_nonzero_sl


import numpy as np

ch = np.loadtxt('stat_0_5778.txt')

show(ch,1300)
ch=select_chi(ch)
ch1= select_nonzero_sl(ch)


show(ch1,1300)

import matplotlib.pyplot as plt
plt.plot(ch1[:,0],ch1[:,6]/2,label='w')


plt.plot(ch1[:,0],ch1[:,1],label='solid')
plt.plot(ch1[:,0],ch1[:,2],label='liquid')
plt.plot(ch1[:,0],ch1[:,3],label='interface')

import matplotlib.pyplot as plt
plt.plot(ch[:,0],ch[:,6]/2,label='w')
