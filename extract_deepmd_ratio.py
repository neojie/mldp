#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
merge to 
@author: jiedeng
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--train_test_ratio","-ttr",type = float,default=3, help="input path file")
parser.add_argument("--OUTCAR","-o",type = str, default = 'OUTCAR', help="OUTCAR name")
parser.add_argument("--deepmd","-d",type = str, default = 'deepmd', help="deepmd folder name")
parser.add_argument("--vaspidx","-vid",type = str, help="idx file, vasp idx, idx[0] >= 1")
parser.add_argument("--idx","-id",type = str, help="idx file, idx[0] >= 0")
parser.add_argument('--test',"-t", default=True, action='store_false',help="Default: save test as set.001? ")

args   = parser.parse_args()


import os
from dpdata import LabeledSystem
import glob
import numpy as np
import shutil

cwd    = os.getcwd()
if args.inputpath:
    from shared_functions import load_paths
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    if isinstance(tmp[0],list):
        paths = []
        for i in range(len(tmp[2])):
            if tmp[1][i] >0: # only keep relax > 0
                paths.append(os.path.join(tmp[2][i],'recal'))
    else:
        paths = [os.path.join(path,'recal') for path in tmp]

else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd)
    paths = [cwd]

def extract_sigma_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            outcar_file.close()
            return sigma
        
def build_fparam(path, outcar, deepmd): # deepmd/..
#    kb = 8.617333262e-5
#    incar = os.path.join(subfolders[0],'INCAR')
    sigma = extract_sigma_outcar(outcar)
    sets=glob.glob(deepmd+"/set*")
    for seti in sets:
        energy=np.load(seti+'/energy.npy')
        size = energy.size
        all_te = np.ones(size)*sigma
        np.save( os.path.join(seti,'fparam.npy'), all_te)

def check_deepmd(path,nsw,outcar,deepmd):
    build_deepmd(path,nsw,outcar,deepmd)
    build_fparam(path, outcar, deepmd)


def build_deepmd(path,nsw,outcar,deepmd):
    ls = LabeledSystem(outcar,fmt='outcar')
    """
    sub_ls = ls.sub_system(idx)
    
    """
    if args.idx:
        print("index file provided")
        idx = np.loadtxt(args.idx).astype(int)
#        ls = ls.sub_system(idx)
    elif (not args.idx) and  args.vaspidx:
        print("vasp index file provided")
        vaspidx=np.loadtxt(args.vaspidx)
        fp  = open(outcar)
        fp.readline()
        nsw_sel = fp.readline()
        if 'nsw_sel' in nsw_sel:
            print('file generated by merge_out.py')
        #    print(nsw_sel)
            tmp     = nsw_sel.split('=')[1].strip().split(' ')
            nsw_sel = [int(tmp_idx) for tmp_idx in tmp]
            idx = []    
            for i in range(len(nsw_sel)):
                if nsw_sel[i] in vaspidx:
                    idx.append(i)
        else:
            print('OUTCAR file generated by VASP')
            idx = vaspidx - 1
    else:
        print("split train and test by ratio {0} : {1}".format(args.train_test_ratio,1))
        train_size = round(len(ls)*(args.train_test_ratio)/(args.train_test_ratio+1))
#        test_size = round(len(ls)*1/(args.train_test_ratio+1))
        idx = np.random.choice(range(len(ls)), train_size, replace=False)
        idx.sort()
    
    idx2 = [i for i in range(len(ls)) if i not in idx] # test
    ls2  = ls.sub_system(idx2) # test
    ls   = ls.sub_system(idx)
        
    deepmd = os.path.join(path,deepmd)
    
    ls.to_deepmd_npy(deepmd,set_size=1000000) # give a *large* value, default is 5000
    if args.test:
        ls2.to_deepmd_npy('test_tmp',set_size=1000000)
        shutil.copytree('test_tmp/set.000',os.path.join(deepmd,'set.001'))
        shutil.rmtree('test_tmp')

def extract_nsw_outcar(outcar):

    """
    length of outcar
    
    """
    outcar_file = open(outcar)
    nsw = 0
    for _ in range(100000000):
        line = outcar_file.readline()
        if 'nsw_tot' in line:
            nsw = float(line.split('=')[1].split()[0])
            outcar_file.close()
            return nsw
        elif 'free  energy   TOTEN' in line:
            nsw += 1
    return nsw

for path in paths:
    print(path)
    if os.path.exists(os.path.join(path,args.deepmd)):
        print("deepmd foler already exist=> skip")
    elif not os.path.exists(os.path.join(path,args.OUTCAR)):
        print("OUTCAR folder do not exists")
    else: # has required OUTCAR folder but no deepmd 
        print("Build {0}".format(args.deepmd))
        nsw =extract_nsw_outcar(os.path.join(path,args.OUTCAR))
        check_deepmd(path,nsw,os.path.join(path,args.OUTCAR), os.path.join(path,args.deepmd))

