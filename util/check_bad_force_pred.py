#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:21:25 2021

@author: jiedeng
"""

import numpy as np
bad_force_threshold = 5
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--detail_file", type=str, 
                        help="The file containing details of energy force and virial accuracy")
parser.add_argument("-bad_force_prediction", "--bf", type=float,default=4, 
                        help="bad force prediction criterion, as difference between force and predictio")   
args   = parser.parse_args()

deepmd_path = '.'
e_test_file = os.path.join(deepmd_path, args.detail_file+".e.out")
f_test_file = os.path.join(deepmd_path, args.detail_file+".f.out")
v_test_file = os.path.join(deepmd_path, args.detail_file+".v.out")
    
f_test = np.loadtxt(f_test_file)
dif = np.abs(f_test[:,3:]-f_test[:,:3])
natoms = 160
print(np.where(dif.max(axis=1)>args.bad_force_prediction)[0]//natoms)
