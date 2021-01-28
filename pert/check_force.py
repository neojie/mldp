#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020
shared function, select what to drop
1-analysis deepmd only, analysis outcar only
outcar only, if no model provided


only consider 
@author: jiedeng
"""


import argparse
import os
import numpy as np
parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file. support txt file, json file, default is cwd")


args   = parser.parse_args()
cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    if len(inputpath) >4 and inputpath[-4:] == 'json':
        print("input is json file, load json")
        import json
        with open(inputpath) as f:
            tmp = json.load(f)
            paths = tmp['training']['systems']
    else:
        from shared_functions import load_paths
        paths = load_paths(inputpath,level='recal')
#    paths = [os.path.join(path,args.suffix) for path in tmp]
else:
    print("No folders point are provided. Use current working path")
#    inputpath = os.path.join(cwd)
    paths = [cwd]

#

for path in paths:
    print('--'*40)
    print(path)

    force=np.load(os.path.join(path,'set.000/force.npy'))
    force1=np.load(os.path.join(path,'set.001/force.npy'))
    print('set.000 and set.001')
    print(min(force.min(axis=1)), max(force.max(axis=1)))
    print(min(force1.min(axis=1)), max(force1.max(axis=1)))
    print('--'*40)