#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 13:54:35 2020

@author: jiedeng
"""

import numpy as np
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--idxs","-ids",type = str, help="idx files, separated by -")
parser.add_argument("--out","-o",type = str, help="idx file, idx[0] >= 0")
args   = parser.parse_args()

paths=args.idxs.split('-')
idx = np.array([])
for path in paths:
    print('start to merge {0}'.format(path))
    tmp=np.loadtxt(path).astype(int)
    print('{0} has {1} indexs'.format(len(tmp)))
    idx = np.union1d(idx,tmp)

print('final file {0} has {1} indexs'.format(args.out,len(idx)))
np.savetxt(args.out,idx,fmt='%d',header=args.idxs)