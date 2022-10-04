#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 20:21:49 2022

@author: jiedeng
"""

import os
import glob
dirs=os.listdir('.')
train = []
test = []
import numpy as np

sets = ['set.000','set.001', 'set.002', 'set.003', 'set.004', \
'set.005', 'set.006', 'set.007', 'set.008', 'set.009', \
'set.010', 'set.011', 'set.012', 'set.013', 'set.014']

for diri in dirs:
    deepmd=os.listdir(diri)[0]
    setf = glob.glob(diri+'/'+deepmd+'/set*')
    
    count = 0 
    for seti in sets:
        seti_path = diri+'/'+deepmd+'/'+seti+'/energy.npy'
        if os.path.exists(seti_path):
            count += len(np.load(seti_path)) 
            
    if len(setf)>1:
        seti= setf[0][:-1]+str(len(setf)-1)
        tmp = len(np.load(seti+'/energy.npy'))
    else:
        tmp = 0
    test.append(tmp)
    train.append(count - tmp)


