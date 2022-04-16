#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:21:25 2021

example
python ~/script/mldp/util/check_efv.py -ct -ip input_v2_compat.json -d m12v3 -bf 20
python check_efv.py -ip input_v2_compat.json -d m12v3 -v -100 300 -vb -ct
@author: jiedeng
"""

import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file. support txt file, json file, default is cwd")
parser.add_argument("--verbose", "-vb", default=True,action='store_false',
                        help="verbose? True show details, False only erroneous ones")

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
            try:
                paths = tmp['training']['systems']
            except:
                paths = tmp['training']['training_data']['systems']
    else:
        from shared_functions import load_paths
        paths = load_paths(inputpath,level='recal')
#    paths = [os.path.join(path,args.suffix) for path in tmp]
else:
    print("No folders point are provided. Use current working path")
#    inputpath = os.path.join(cwd)
    paths = [cwd]

count = 0
for path in paths:
    # if args.verbose:
    #     print('--'*40)
    #     print(path)
    #     print('--'*40)
    
    if os.path.exists(os.path.join(path,'set.001')):
        pass
    else:
        print('--'*40)
        print(path)
        print('--'*40)
        print(np.load(os.path.join(path,'set.000/energy.npy')).shape[0])