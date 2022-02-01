#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:21:25 2021

@author: jiedeng
"""

import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file. support txt file, json file, default is cwd")

parser.add_argument("-d", "--detail_file", type=str, 
                        help="The file containing details of energy force and virial accuracy")
parser.add_argument("--bad_force_prediction", "-bf", type=float, 
                        help="Default: None, recommended 4, bad force prediction criterion, as difference between force and predictio")   
parser.add_argument("--force_range",'-f', type=float,nargs="+", 
                        help="Default: None, force (eV/A) range") 
parser.add_argument("--energy_peratom_range",'-e',type=float,nargs="+", 
                        help="Default: None, energy (eV/atom) range")   
parser.add_argument("--virial_range",'-v', type=float,nargs="+", 
                        help="Default: None, virial (GPa) range peratom")  
parser.add_argument('--check_test',"-ct", default=True, action='store_false',
                    help="Default: check test, if false, check train")

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
    print('--'*40)
    print(path)
    print('--'*40)
    ## check the source data, some only has deepmd, not deepmd_relax2
    ## and vice versa. deepmd_relax2 is for relax process only
    ## to do the statisitics, we only care where there is deepmd
    try:
        assert os.path.exists(os.path.join(path,args.deepmd))  # deepmd path must exist
        deepmd_path = os.path.join(path,args.deepmd)
    except:
        assert os.path.exists(path) # this gives the fixability 
        deepmd_path = path
        
    
    eV_A3_2_GPa  = 160.21766208 # 1 eV/Å3 = 160.2176621 GPa
    natoms = np.loadtxt(os.path.join(deepmd_path,'type.raw')).shape[0]  #assume natoms does not change
    import glob
    sets=glob.glob(os.path.join(deepmd_path,'set.*'))
    boxs = np.load(os.path.join(sets[0],'box.npy')) #box may shape may vary from frame to frame => NPT, MPMT
    box  = boxs[0]
    a    = box[:3]
    b    = box[3:6]
    c    = box[6:]
    vol  = np.dot(c,np.cross(a,b))
        
    natoms = np.loadtxt(os.path.join(deepmd_path,'type.raw')).shape[0]  #assume natoms does not change
    
    if args.check_test:
        e_test_file = os.path.join(deepmd_path, args.detail_file+".e.out")
        f_test_file = os.path.join(deepmd_path, args.detail_file+".f.out")
        v_test_file = os.path.join(deepmd_path, args.detail_file+".v.out")
    else:
        e_test_file = os.path.join(deepmd_path, args.detail_file+".e.tr.out")
        f_test_file = os.path.join(deepmd_path, args.detail_file+".f.tr.out")
        v_test_file = os.path.join(deepmd_path, args.detail_file+".v.tr.out")        
    
    e_test = np.loadtxt(e_test_file)
    f_test = np.loadtxt(f_test_file)
    v_test = np.loadtxt(v_test_file)
    v_gpa_test = v_test/vol*eV_A3_2_GPa 
    
#    print(v_gpa_test)
    
    bad_force_exclude_idx = []
    energy_exclude_idx = []
    force_exclude_idx  = []
    virial_exclude_idx = []
    
    out = []
    
    if args.bad_force_prediction:
        dif    = np.abs(f_test[:,3:]-f_test[:,:3])
        bad_force_exclude_idx = np.where(dif.max(axis=1)>args.bad_force_prediction)[0]//natoms
        out=np.union1d(out,bad_force_exclude_idx)
    
        
    if args.energy_peratom_range:
        max_e = max(args.energy_peratom_range)*natoms
        min_e = min(args.energy_peratom_range)*natoms
        tmp1= np.where(e_test[:,1]>max_e)[0]
        tmp2= np.where(e_test[:,1]<min_e)[0]
        energy_exclude_idx  = np.union1d(tmp1,tmp2)
        out=np.union1d(out,energy_exclude_idx)
    
    if args.force_range:
        max_f = max(args.force_range)
        min_f = min(args.force_range)
        tmp1= np.where(f_test[:,:3].max(axis=1)>max_f)[0]
        tmp2= np.where(f_test[:,:3].min(axis=1)<min_f)[0]
        force_exclude_idx  = np.union1d(tmp1,tmp2)//natoms
        out=np.union1d(out,force_exclude_idx)
    
    if args.virial_range:
        max_v = max(args.virial_range)
        min_v = min(args.virial_range)
        tmp1= np.where(v_gpa_test[:,:9].max(axis=1)>max_v)[0]
        tmp2= np.where(v_gpa_test[:,:9].min(axis=1)<min_v)[0]
        virial_exclude_idx  = np.union1d(tmp1,tmp2)
        out=np.union1d(out,virial_exclude_idx)
    if len(out)>0:
        print("Frames dis-satfitisfy", out.astype(int)) 
    else:
        print("**good**")



