#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 16:35:11 2020

generate poscar from lammps files
@author: jiedeng
"""
import numpy as np
import os
from shutil import copy
import re
from shared_functions import load_paths, check_outcar_done
import dpdata
import glob

def lmp2pos(path,step,copybs=False):
    """
    build POSCAR based on XDATCAR given in the path
    copy INCAR, POTCAR, KPOINTS into the same folder
    
    """
    lmps1=glob.glob(os.path.join(path,'*lammps*'))
    lmps2=glob.glob(os.path.join(path,'*dump*'))
    if lmps2 is []:
        print("No dump file found")
    elif not(lmps2 is []):
        lmp = lmps2[0]    
    elif not(lmps1 is []):
        lmp = lmps1[0]
    else:
        print("No dump file and lammps file found")
        
    path = os.path.join(path,'recal')

    ls=dpdata.System(lmp,fmt='lammps/dump')
    for i in range(0,len(ls),step):
        try:
            os.mkdir(os.path.join(path,str(i+1))) # return None
        except:
            print("Folder {0} already exists,skip making".format(i))
        target_path    = os.path.join(path,str(i+1))
        ls_target_path = os.listdir(target_path)
                    
        if 'OUTCAR'in ls_target_path and check_outcar_done(os.path.join(target_path,'OUTCAR')):
            print('calculation done')
        elif 'INCAR'   in ls_target_path and \
             'POTCAR'  in ls_target_path and \
             'POSCAR'  in ls_target_path and \
             'KPOINTS' in ls_target_path:
            pass
        else:
            ls.to_vasp_poscar(os.path.join(target_path,'POSCAR'),frame_idx=i)  
            copy(os.path.join(inputfile,'INCAR'),target_path)
            copy(os.path.join(inputfile, 'KPOINTS'),target_path)
            copy(os.path.join(inputfile, 'POTCAR'),target_path)  
        run(cwd,target_path)

from subprocess import call
def run(cwd,target_path):
    os.chdir(target_path)
    sub_file = os.path.join(inputfile,'sub_vasp.sh')
    call("qsub {0}".format(sub_file), shell=True)
#    call("sbatch sub_script.sh", shell=True)
    os.chdir(cwd)
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--inputfile","-if",help="input files for vasp cal")
parser.add_argument("--step","-s",default=1,type=int,help="step")
parser.add_argument("--range","-r",default='1_-1',type=str,help="step")

args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath)
else:
    print("input path are provided.")
    paths = [cwd]

if args.inputfile:
    print("Check files in {0}  ".format(args.inputfile))
    inputfile = args.inputfile
else:
    print("No folders point are provided. Use default value folders")
    inputfile = os.path.join(cwd,'inputs')
    
for path in paths:
    print('###',path)

    if os.path.exists(os.path.join(path,'recal')):
        print('*** SKIP',path,'has recal folder')
    else:
        print('--> Build recal',path)
        os.mkdir(os.path.join(path,'recal'))      
        ### create jobs
        lmp2pos(path,args.step,copybs = True)
#        dsq_jobs(path,sel_nsw)
#        run_dsq(cwd,path)
    

