#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:53:36 2021

select model stat results

@author: jiedeng
"""

import numpy as np
import argparse
import os
parser = argparse.ArgumentParser()

parser.add_argument("-d", "--detail_file", type=str, 
                        help="The file containing details of energy force and virial accuracy")
parser.add_argument('--exclude',"-e", type=int,nargs="+",help="manually exclude indexs") 
args   = parser.parse_args()


exclude = args.exclude
natoms = 160
deepmd_path = '.'
e_test_file = os.path.join(deepmd_path, args.detail_file+".e.out")
f_test_file = os.path.join(deepmd_path, args.detail_file+".f.out")
v_test_file = os.path.join(deepmd_path, args.detail_file+".v.out")

e_test = np.loadtxt(e_test_file)
f_test = np.loadtxt(f_test_file)
v_test = np.loadtxt(v_test_file)

idx_new = [i for i in range(len(e_test)) if i not in exclude]
exclude_f = []
for i in exclude:
    for j in range(natoms):
        exclude_f.append(i*natoms +j)

idx_new_f = [i for i in range(len(e_test)*natoms) if i not in exclude_f]

e_sel = e_test[idx_new,:]
f_sel = f_test[idx_new_f,:]
v_sel = v_test[idx_new,:]

import os
os.rename(e_test_file,e_test_file+'_org')
os.rename(f_test_file,f_test_file+'_org')
os.rename(v_test_file,v_test_file+'_org')

np.savetxt(e_test_file,e_sel)
np.savetxt(f_test_file,f_sel)
np.savetxt(v_test_file,v_sel)
