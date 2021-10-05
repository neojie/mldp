#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 13:30:00 2020

@author: jiedeng
"""

### select frames from test test
import os  
import numpy as np
import argparse

parser = argparse.ArgumentParser()
#l1=dpdata.System('/Users/jiedeng/Documents/ml/deepmd-kit/my_example/tail_3.dump')
parser.add_argument("--cutoff","-c",type=float,default=4,help="cutoff of force difference")
args   = parser.parse_args()

f_test_file = 'mm10'+".f.out"
f_test = np.loadtxt(f_test_file)
dif_f_test = np.abs(f_test[:,:3] - f_test[:,3:])
cutoff = args.cutoff
anomalous_idx = np.unique(np.where(dif_f_test>cutoff)[0]//160)
npys = ['energy.npy', 'box.npy',  'coord.npy' , 'force.npy' , 'fparam.npy' , 'virial.npy']

energy = np.load('set.001/energy.npy')
idxs=[i for i in range(len(energy)) if i not in anomalous_idx]

os.mkdir('reduced.001')
for npy in npys:
    tmp = np.load('set.001/'+npy)[idxs]
    np.save('reduced.001/'+npy,tmp)

print('{0} out of {1} frame are DEselected'.format(len(anomalous_idx), len(energy)))
np.savetxt('reduced.001/anomalous_idx.txt',anomalous_idx,fmt='%d',header='np.unique(np.where(dif_f_test>{0})[0]//160)'.format(cutoff))
