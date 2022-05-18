#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 21:20:52 2022

#  id    H           O           Mg           Si           chi
0    0 212 4    414 520 254    133 132 95    139 144 77    0.0094

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
parser.add_argument("--num_interface_w","-nw",type=int,help="Default: larger, btetter, must set!")
parser.add_argument("--mode","-m",help="Default: mass mode, k or mass, must set!")


# parser.add_argument("--elements","-ele",nargs="+",help="elements to be analyzed") xyz file has that

args   = parser.parse_args()


from gds_analyzer import GDSAnalyzer


import sys,os
try:
    sys.path.insert(1, '/Users/jiedeng/GD/papers/pv7_h/partition/codes/')
except:
    print("run in local")
try:
    sys.path.insert(1, '/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/mldp/similiarity/')
except:
    print("run in local")

if args.file:
    xyz = args.file
else:
    cwd = os.path.abspath(os.curdir)
    xyz     = os.path.join(cwd,'merge.xyz')
    
ana=GDSAnalyzer(xyz = xyz, begin = args.begin, end = args.end, mode = args.mode, project_axis= args.project_axis, step = args.step,nw=args.num_interface_w)

