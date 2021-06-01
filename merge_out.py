#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
merge to 
@author: jiedeng
"""
from shared_functions import load_paths, merge_sel, remove_recal_traj_files
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
#parser.add_argument("--mode","-m",default = 'sel',type=str,help="input path file")
parser.add_argument("--outcar","-o",type=str,help="name of outcar")
parser.add_argument("--remove","-r",type=str,help="if input, remove subfolders") # intentionally do not input default value

args   = parser.parse_args()
cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    if isinstance(tmp[0],list): # relax [[x],[x],[x]], regular [x,x,x]
        paths = []; relax_steps = tmp[1]; relax_paths = tmp[2]
        for i in range(len(relax_steps)):
            if relax_steps[i]>0:
                paths.append(os.path.join(relax_paths[i],'recal'))
    else:
        paths = [os.path.join(path,'recal') for path in tmp]
else:
    print("No folders point are provided. Use default value folders")
    paths = [cwd]
    

for i in range(len(paths)):
    path = paths[i]
    print(path)
    merge_sel(path, args.outcar)
    if os.path.exists(os.path.join(path,args.outcar)) and args.remove:
        remove_recal_traj_files(path)
            
    
            

