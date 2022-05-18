#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:10:33 2022
submit a batch of analyze scripts
@author: jiedeng
"""
import argparse
parser = argparse.ArgumentParser()

## must set
parser.add_argument("--num_interface_w","-nw",type=int,help="Default: larger, btetter, must set!")
parser.add_argument("--mode","-m",help="Default: mass mode, must set!")

# be alert the possible change, for instance, ppv
parser.add_argument("--project_axis","-p",default=2,type=int,help="default 2(z), 0,1,2 => x,y,z")

# you may want to adjust
parser.add_argument("--interval","-i",default=30,type=int,help="interval")
parser.add_argument("--file","-f",type=str,help="path to xyz file to analyze, default is merge.xyz in the cwd")

# normally leave them alone
parser.add_argument("--end","-e",type=int,help="end index, if file length is 100, then 100, not 99")
parser.add_argument("--begin","-b",default=0,type=int,help="begining index, default: 0 ")
parser.add_argument("--run","-r",default=True,action='store_false',help="submit? default : Yes")
parser.add_argument("--step","-s",default=1,type=int,help="step")


args   = parser.parse_args()

# import numpy as np
from subprocess import call


interval  = args.interval
beg = args.begin



print('--'*20)
print('summary of setting')
print('interface width = ',args.num_interface_w)
print('--'*20)
string = """#!/bin/bash
#$ -cwd
#$ -o out.$JOB_ID
#$ -j y
#$ -pe shared 1
#$ -l h_rt=24:00:00
. "/u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh"
conda activate /u/home/j/jd848/project-lstixrud/dpkit2;
. /u/local/Modules/default/init/modules.sh
module load gcc/8.3.0;module load intel/2020.4;module load cmake

export PLUMED_KERNEL=/u/home/j/jd848/project-lstixrud/plumed/lib/libplumedKernel.so
export PATH=/u/home/j/jd848/project-lstixrud/plumed/bin:$PATH
export LD_LIBRARY_PATH=/u/home/j/jd848/project-lstixrud/plumed/lib:$LD_LIBRARY_PATH
"""
import os

if args.file:
    xyz = args.file
else:
    cwd = os.path.abspath(os.curdir)
    xyz     = os.path.join(cwd,'merge.xyz')

if args.end:
    end = args.end
else:
    import MDAnalysis as mda
    mda_xyz = mda.Universe(xyz)
    end = len(mda_xyz.trajectory)
    
# if os.path.exists('log.sub'):
#     log = open('log.sub','a')
# else:
    # log = open('log.sub','w')
log = open('log.sub','w')
for i in range(beg,end,interval):
    file = open('sub_{0}'.format(i//interval),'w')
    file.writelines(string)
    if i+interval < end:
        endidx = i+interval
    else:
        endidx = end
    file.writelines('python ~/script/mldp/similarity/stat.py -sh -b {0} -e {1} -p {2} -f {3} -m {4} -s {5} -nw {6}'.format(i,endidx,args.project_axis, xyz,args.mode,args.step, args.num_interface_w))
    log.writelines('stat_{0}_{1}.txt\n'.format(i,endidx-1))
    file.close()
    if args.run:
        call("/u/systems/UGE8.6.4/bin/lx-amd64/qsub sub_{0}".format(i//interval),shell=True)        
log.close()
