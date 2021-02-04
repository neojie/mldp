#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:21:25 2021

@author: jiedeng
"""

import numpy as np
f_test = np.loadtxt('om8.f.out')
dif = np.abs(f_test[:,3:]-f_test[:,:3])
natoms = 160
print(np.where(dif.max(axis=1)>5)[0]//natoms)