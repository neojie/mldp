#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 19:35:28 2021

@author: jiedeng
"""

import glob
import os
from shutil import copy
from subprocess import call
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--inputfile","-if",help="input files for vasp cal")
#parser.add_argument("--step","-s",default=1,type=int,help="step")
#parser.add_argument("--range","-r",type=str,help="0-2, means from 0 to 2, default is for all folders")
parser.add_argument("--run_vasp","-rv",default=True,action='store_false',help="run vasp?, default without input is Yes")

args   = parser.parse_args()

def run(target_path):
    os.chdir(target_path)
    sub_file = os.path.join(args.inputfile,'sub_vasp.sh')
    call("/u/systems/UGE8.6.4/bin/lx-amd64/qsub {0}".format(sub_file), shell=True)
#    call("bash {0}".format(sub_file), shell=True)
    os.chdir(cwd)
    
pos = glob.glob('p*.vasp')
cwd = os.getcwd()
    
for i in range(len(pos)):
    os.mkdir(str(i))
    copy('p'+str(i)+'.vasp',str(i)+'/POSCAR')
    copy(os.path.join(args.inputfile,'INCAR'),str(i))
    copy(os.path.join(args.inputfile, 'KPOINTS'),str(i))
    copy(os.path.join(args.inputfile, 'POTCAR'),str(i))  
    if args.run_vasp:        
        run(str(i))
    

        

    