#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 22:41:45 2020

@author: jiedeng
"""

import numpy as np
import os
import glob
import argparse
parser = argparse.ArgumentParser()
#l1=dpdata.System('/Users/jiedeng/Documents/ml/deepmd-kit/my_example/tail_3.dump')
parser.add_argument("--sets","-s",type=str,nargs="+",help="set files, e.g., ,x.001 x.002")
parser.add_argument("--os","-os",type=str,default='merge.001',help="merged file, default: merge.001")

args   = parser.parse_args()

sets = args.sets
count = 0
#npys = ['box.npy',  'coord.npy' , 'energy.npy' , 'force.npy' , 'fparam.npy' , 'virial.npy']
npys =glob.glob(sets[0]+'/*npy')


out = []
for npy in npys:
    out.append(np.load(npy))

try:
    for i in range(1,len(sets)):
        SET = sets[i]
        npys = glob.glob(SET+'/*npy')
        for j in range(len(npys)):
            out[j] = np.concatenate((out[j],np.load(npys[j])),axis=0)
except:
    print("sets must >1")

set_name = npys[0].split('/')[0]
tar_name = args.os


if not os.path.exists(tar_name):
    os.mkdir(tar_name)

for npy, ou  in zip(npys,out):
    name = npy.replace(set_name,tar_name)
    np.save(name,ou)