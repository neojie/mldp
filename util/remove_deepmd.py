#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 21:05:29 2020

@author: jiedeng
"""
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
import shutil
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
    
def build_deepmd(path,nsw):
    ls = LabeledSystem(os.path.join(path, 'outcar'),fmt='outcar')
    deepmd = os.path.join(path,'deepmd')
    if nsw <= 4: # we know nsw must > 100
        set_size = 1
        print("{0} has only {1}".format(path,nsw))
    if nsw > 4:
        set_size = nsw//4  # 25% used as , but if say 82, then 20, 20, 20, 2, too less
    ls.to_deepmd_npy(deepmd,set_size=set_size)
    
for path in paths:
    print(path)
    if os.path.exists(os.path.join(path,'post_recal_done')):
        if os.path.exists(os.path.join(path,'deepmd')):
            shutil.rmtree(os.path.join(path,'deepmd'))
        if os.path.exists(os.path.join(path,'outcar')):
            os.remove(os.path.join(path,'outcar'))
            print('remove done')
    else:
        print('No post_recal done?')


