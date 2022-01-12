#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 16:35:11 2020

generate poscar from deepmd files

if recal exist, check if the NSW is selected?
@author: jiedeng
"""
import os
from shutil import copy
#from shared_functions import load_paths
import dpdata
#import glob

def lmp2pos(ls,sel_nsw,copybs=False):
    """
    build POSCAR based on deepmd given in the path
    copy INCAR, POTCAR, KPOINTS into the same folder
    
    """
#    lmps1=glob.glob(os.path.join(path,'*lammps*'))
#    lmps2=glob.glob(os.path.join(path,'*dump*'))
#    if lmps2 is []:
#        print("No dump file found")
#    elif not(lmps2 is []):
#        lmp = lmps2[0]    
#    elif not(lmps1 is []):
#        lmp = lmps1[0]
#    else:
#        print("No dump file and lammps file found")
        
#    path = os.path.join(path,'recal')

#    ls=dpdata.System(lmp,fmt='lammps/dump')
#    recal_path = os
    if  sel_nsw is None:
        sel_nsw = range(0,len(ls),args.step)       
    else:
        sel_nsw = sel_nsw
    for i in sel_nsw:
        print(i)
        if os.path.exists(os.path.join(recal_path,str(i+1))):
            print("Folder {0} already exists,skip making".format(i))
        else:
            os.mkdir(os.path.join(recal_path,str(i+1))) # return None
            target_path    = os.path.join(recal_path,str(i+1))                   
            ls.to_vasp_poscar(os.path.join(target_path,'POSCAR'),frame_idx=i)  
            copy(os.path.join(inputfile,'INCAR'),target_path)
            copy(os.path.join(inputfile, 'KPOINTS'),target_path)
            copy(os.path.join(inputfile, 'POTCAR'),target_path)  
            if run_vasp:
                print("run vasp",target_path)
                run(cwd,target_path)

from subprocess import call
def run(cwd,target_path):
    os.chdir(target_path)
    sub_file = os.path.join(inputfile,'sub_vasp.sh')
    call("/u/systems/UGE8.6.4/bin/lx-amd64/qsub {0}".format(sub_file), shell=True)
#    call("bash {0}".format(sub_file), shell=True)
    os.chdir(cwd)
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--deepmd","-d",help="deepmd path file")
parser.add_argument("--inputfile","-if",help="input files for vasp cal, default is cwd+inputs, if input, please input the absolute path")
parser.add_argument("--step","-s",default=1,type=int,help="step")
parser.add_argument("--range","-r",type=str,help="0-2, means from 0 to 2, default is for all folders")
parser.add_argument("--run_vasp","-rv",help="run vasp?, default without input is Yes")

args   = parser.parse_args()
print(args.run_vasp)
cwd    = os.getcwd()
if args.deepmd:
    print("Check files in {0}  ".format(args.deepmd))
#    inputpath = args.inputpath
#    paths = load_paths(inputpath)
    ls = dpdata.System(args.deepmd,fmt='deepmd/npy')
else:
    print("deepmd path is not provided .")
    ls = dpdata.System('.',fmt='deepmd/npy')

if args.inputfile:
    print("Check files in {0}  ".format(args.inputfile))
    inputfile = args.inputfile
else:
    print("No folders point are provided. Use default value folders")
    inputfile = os.path.join(cwd,'inputs')

# change inputfile to absolute path, chdir occurs when submitting job

inputfile = os.path.abspath(inputfile)
sel_nsw = None
if args.range:
    tmp     = [int(i) for i in args.range.split('-')]
    sel_nsw = range(tmp[0],tmp[1],args.step)
run_vasp = True
if args.run_vasp:
    run_vasp = False
   
#for path in paths:
#    print('###',path)
recal_path = os.path.join(os.getcwd(),'recal')
try:
    os.mkdir(recal_path)     
except:
    print('***recal exists in',recal_path)

lmp2pos(ls,sel_nsw,copybs = True)
    

