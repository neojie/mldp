#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 13:30:00 2020

@author: jiedeng
"""

### select frames from test test
import os  
import numpy as np
f_test_file = 'mm10'+".f.out"
f_test = np.loadtxt(f_test_file)
dif_f_test = np.abs(f_test[:,:3] - f_test[:,3:])
anomalous_idx = np.unique(np.where(dif_f_test>4)[0]//160)
npys = ['energy.npy', 'box.npy',  'coord.npy' , 'force.npy' , 'fparam.npy' , 'virial.npy']

energy = np.load('set.001/energy.npy')
idxs=[i for i in range(len(energy)) if i not in anomalous_idx]

os.mkdir('reduced.001')
for npy in npys:
    tmp = np.load('set.001/'+npy)[idxs]
    np.save('reduced.001/'+npy,tmp)
np.savetxt('reduced.001/anomalous_idx.txt',anomalous_idx,fmt='%d',header='np.unique(np.where(dif_f_test>4)[0]//160)')
