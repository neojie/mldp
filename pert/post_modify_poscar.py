#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 20:23:32 2021
use 
merge_out.py
extract_deepmd

1- dp test with model prefix
2- summarize dp test results compare energy, force, and virial

example:
    python ~/script/mldp/check_nbands_nelm.py -ip all
    python ~/script/mldp/merge_out.py -o OUTCAR -r y
    python ~/script/mldp/extract_deepmd.py -bs 1000 -o OUTCAR -d deepmd -t 
(deepmd-kit-gpu)[jd848@n7214 deepmd]$ dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/mpmt/mm12/pv_Jan8_cpu.pb -d mm12

dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/mpmt/mm12/pv_Jan8_cpu.pb -d mm12
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/extreme_filtered/re4/pv-cpu.pb -d extreme4  
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/extreme_filtered/re7/pv_gpu.pb -d extreme7   
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/pert/mp1/pv.pb -d mp1 
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/onlymelt/om2/pv.pb -d om2


mgsi done
@author: jiedeng
"""


import argparse
parser = argparse.ArgumentParser("Test model on poscars where atoms are manually palced and interatomic distances are controlled ")
parser.add_argument("--index","-i",nargs="+", help="index of atom of interest, if not provided, read from log file")
parser.add_argument("--log","-log", help="log file")
parser.add_argument("--test_folder","-tf",default='deepmd', help="where is dp test located?")
parser.add_argument("--model_prefix","-mp", help="model prefix")
parser.add_argument("--natoms","-na",type=int,default = 20, help="log file")
args = parser.parse_args()

import numpy as np
import os

if args.index:
    atom1, atom2 = args.index[0], args.index[1]
else:
    if args.log:
        tmp = np.loadtxt(args.log).astype(int)
    else:
        tmp = np.loadtxt('log').astype(int)
    atom1, atom2 = tmp[0], tmp[1]
    
prefix = args.model_prefix
test_folder = args.test_folder
e_dir = os.path.join(test_folder,prefix+'.e.out')
f_dir = os.path.join(test_folder,prefix+'.f.out')
v_dir = os.path.join(test_folder,prefix+'.v.out')


## calc distance
coord_dir = os.path.join(os.path.join(test_folder,'set.000'),'coord.npy')
box_dir   = os.path.join(os.path.join(test_folder,'set.000'),'box.npy')

coord     = np.load(coord_dir)
box = np.load(box_dir)[0]
box = np.array([box[0], box[4], box[8]])
coord_atom1 = coord[:,(atom1*3):(atom1*3+3)]
coord_atom2 = coord[:,(atom2*3):(atom2*3+3)]

#cal distance 
diffs = coord_atom2 - coord_atom1
cond = np.where(diffs > box / 2.)
diffs[cond] -= box[cond[1]]
cond = np.where(diffs < -box / 2.)
diffs[cond] += box[cond[1]]
sqdist = np.square(diffs).sum(axis=1)
dist = np.sqrt(sqdist)     


e=np.loadtxt(e_dir)
f=np.loadtxt(f_dir)
v=np.loadtxt(v_dir)

# len(e), num of frames
f_idx = np.array(range(len(e)))*args.natoms + atom2
f2 = f[f_idx,:] #force felt by atom2
f2_vasp = f2[:,:3]; mag_f2_vasp = np.sqrt(np.sum(f2_vasp**2,axis=1))
f2_nn   = f2[:,3:]; mag_f2_nn = np.sqrt(np.sum(f2_nn**2,axis=1))

import matplotlib.pyplot as plt
fig,ax = plt.subplots(2,1,figsize=(4,6),sharex=True)
ax[0].plot(dist,mag_f2_nn,'o',label='nn')
ax[0].plot(dist,mag_f2_vasp,'+',label='vasp')

ax[1].plot(dist,e[:,1],'o',label='nn')
ax[1].plot(dist,e[:,0],'+',label='vasp')

ax[0].set_ylabel('force (eV/A)')
ax[1].set_ylabel('energy (eV)')
ax[1].set_xlabel('interatomic distance (A)')
ax[0].grid()
ax[1].grid()
ax[0].legend()
ax[1].legend()
plt.minorticks_on()
plt.show()






