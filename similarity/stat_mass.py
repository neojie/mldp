#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 21:20:52 2022

@author: jiedeng
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--begin","-b",default=0,type=int,help="begining index, default: 0 ")
parser.add_argument("--end","-e",type=int,help="end index, default is end of file")
parser.add_argument("--project_axis","-p",default=2,type=int,help="default 2(z), 0,1,2 => x,y,z")
parser.add_argument("--step","-s",default=1,type=int,help="step")
parser.add_argument("--file","-f",type=str,help="path to xyz file to analyze, default is merge.xyz in the cwd")
parser.add_argument("--show","-sh",default=True,action='store_false',help="Default: show results")
args   = parser.parse_args()


import MDAnalysis as mda
import sys,os
try:
    sys.path.insert(1, '/Users/jiedeng/GD/papers/pv7_h/partition/codes/')
except:
    print("run in local")
try:
    sys.path.insert(1, '/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/mldp/similarity/')
except:
    print("run in local")

from stat_lib_mass import analyze, show
import ase.io
if args.file:
    xyz = args.file
else:
    cwd = os.path.abspath(os.curdir)
    xyz     = os.path.join(cwd,'merge.xyz')

project_axis = args.project_axis #0,1,2 => x,y,z
start_idx = args.begin # 972   ## omit the first 100 ps, if D is 9e-10 m^2/s, this means 3 A of H diffusion
end_idx    = args.end#2
#----------------------------------------------------------------------#
alpha = 2.5
ase_xyz = ase.io.read(xyz,index=':') 
mda_xyz = mda.Universe(xyz)
# length  =  len(ase_xyz)
if args.end:
    end = args.end
else:
    end = len(ase_xyz)
ch = analyze(args.begin,end, xyz,ase_xyz,mda_xyz, project_axis,alpha, step=args.step, save=True,name='stat_{0}_{1}.txt'.format(args.begin,end-1), assert_chi = False)
if args.show:
    show(ch)
    
    

