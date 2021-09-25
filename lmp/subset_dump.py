#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 21:34:08 2021
build a subset of dump file
@author: jiedeng
"""



import argparse
description="""
        Only the absolute idx is considered. The stride is ignored. If idx is 
        absolute timestep, e.g, 0, 500, 1000. The option of stride can be turned on
        """
print(description)

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type = str, default = 'dump.0', help="dump file name, default: dump.0")
parser.add_argument("--idx","-id",type = str, help="idx file, idx[0] >= 0")
parser.add_argument('--id',"-i", type=int,nargs="+",help="manually input indexs")
parser.add_argument('--stride',"-s", default=False, action='store_true')

args        = parser.parse_args()

import numpy as np
from dump import dump
import copy


d=dump(args.file)
dd=copy.copy(d)
if args.idx == None:
    subset = args.id
else:
    subset = np.loadtxt(args.idx).astype(int)

if args.stride:
    stride = d.snaps[1].time - d.snaps[0].time
    subset = np.array(subset)//stride

snaps =[]
for i in subset:
    snaps.append(dd.snaps[i])
dd.snaps = snaps
dd.nsnaps = dd.nselect = len(subset)
dd.write('subset.dump')

