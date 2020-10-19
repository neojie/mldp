#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
merge to 
@author: jiedeng
"""
from shared_functions import load_paths
from subprocess import call
import os
import argparse
from dpdata import LabeledSystem
import glob
import numpy as np
import re
# check if exist 

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,'recal') for path in tmp]
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd)
    paths = [cwd]
    
def build_outcar(path):
    os.chdir(path)
    if os.path.exist('post_recal_done'):
        call('for i in *;do cat $i/OUTCAR >>outcar;done', shell=True)
    os.chdir(cwd)

def build_fparam(path): # deepmd/..
    kb = 8.617333262e-5
    incar = os.path.join(subfolders[0],'INCAR')
    with open(incar) as incar_file:
        for line in incar_file:
            if 'SIGMA' in line:
               sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            if 'TEBEG' in line:
               tebeg = float(re.sub('[^0-9]', '', line.split('=')[1]))
    if abs(kb*tebeg - sigma)> 0.01:
        print(path,'SIGMA wrong')
        raise ValueError()
    else:
        deepmd = os.path.join(path,'deepmd')
        sets=glob.glob(deepmd+"/set*")
        for seti in sets:
            energy=np.load(seti+'/energy.npy')
            size = energy.size
            all_te = np.ones(size)*sigma
            np.save( os.path.join(seti,'fparam.npy'), all_te)

def check_deepmd(path,nsw):
    build_deepmd(path,nsw)
    build_fparam(path)

def best_size(tot):
    diff = 100000
    for j in range(3,6):## to be 3, 4, 5 parts
        if tot%j == 0 :
            tmp = 0
        else:
            tmp = tot//j - tot%( tot//j )
        print(tmp)
        if tmp < diff:
            diff = tmp
            size = tot//j
    return size,diff

def build_deepmd(path,nsw):
    ls = LabeledSystem(os.path.join(path, 'outcar'),fmt='outcar')
    deepmd = os.path.join(path,'deepmd')
    if nsw <= 4: # we know nsw must > 100
        set_size = 1
        print("{0} has only {1}".format(path,nsw))
    if nsw > 4:
        set_size,_ = best_size(nsw)  # 25% used as , but if say 82, then 20, 20, 20, 2, too less
    ls.to_deepmd_npy(deepmd,set_size=set_size)

for path in paths:
    print(path)
    if os.path.exists(os.path.join(path,'post_recal_done')):
        if not os.path.exists(os.path.join(path,'deepmd')):
            subfolders = [f.path for f in os.scandir(path) if f.is_dir()]
            nsw        = len(subfolders)  # only runs are in folder        
            outcar     = os.path.join(path,'outcar')
            if os.path.exists(outcar):
                print("OUTCAR already exists")
                raise ValueError
            for folder in subfolders:
#                print()
                outcarInfolder = os.path.join(folder,'OUTCAR')

                os.system("cat {0} >> {1}".format(outcarInfolder, outcar))
            check_deepmd(path,nsw)
        else:
            print("deepmd foler already exist")            
    else:
        print("not run post_recal done yet!")

