#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 23:39:38 2020

@author: jiedeng
"""

import glob
import numpy as np
import shutil
import os
import argparse
from shared_functions import load_paths

parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--deepmd","-d",type = str, default = 'deepmd', help="deepmd folder name")
parser.add_argument("--deepmd_scale","-s",type = str, default = 'deepmd_scale', help="scaled deepmd folder name")
parser.add_argument("-S", "--set_prefix", default="set", type=str, help="The set prefix")
parser.add_argument("-or", "--override", type=str, help="The set prefix")
args = parser.parse_args()
cwd  = os.getcwd()

set_prefix = args.set_prefix


if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath,level='recal')
    if isinstance(tmp[0],list):
        paths = tmp[0]
    else:
        paths = tmp

else:
    print("No folders point are provided. Use current working folders")
    paths = [cwd]
    

def Fstd(log_sigma,log_vol,a=0.52707802,b=8.28553082,c=-0.80170509):
    # a=0.52707802,b=8.28553082,c=-0.80170509
    # a=0.53614978,b=8.42219759,c=-0.82133226 used before
    log_fstd = a*log_sigma + c*log_vol + b 
    return np.exp(log_fstd)

def scale_model(vol,sigma):    
    log_vol = np.log(vol)
    log_sigma = np.log(sigma)   
    return Fstd(log_sigma,log_vol)


for path in paths:

    deepmd_path = os.path.join(path,args.deepmd)    
    scaled_path = os.path.join(path,args.deepmd_scale)#deepmd_path.replace(deepmd,scal)
    if args.override and os.path.exists(scaled_path):
        shutil.rmtree(scaled_path)
    shutil.copytree(deepmd_path, scaled_path)   
    dirs        = glob.glob (os.path.join(deepmd_path, set_prefix + ".*"))
    dirs_scaled = glob.glob (os.path.join(scaled_path, set_prefix + ".*"))   
    
    for i in range(len(dirs)):
        dire = dirs[i]
        dire_scaled = dirs_scaled[i]
    
        energy = np.load(os.path.join(dire,'energy.npy'))
        force  = np.load(os.path.join(dire,'force.npy'))
        virial = np.load(os.path.join(dire,'virial.npy'))
        box    = np.load(os.path.join(dire,'box.npy'))
        coor   = np.load(os.path.join(dire,'coord.npy'))
        fparam = np.load(os.path.join(dire,'fparam.npy'))
        vol    = box[0][0]*box[0][4]*box[0][8]   # vol format  [10.626586  0. 0.   0.  10.626586  0.    0. 0.   10.626586]
        scale  = scale_model(vol, fparam[0])
        
        energy_scaled = energy/scale
        force_scaled  = force/scale
        virial_scaled = virial/scale
        
        np.save(os.path.join(dire_scaled,'energy.npy'),energy_scaled)
        np.save(os.path.join(dire_scaled,'force.npy'),force_scaled)
        np.save(os.path.join(dire_scaled,'virial.npy'),virial_scaled)       
    
