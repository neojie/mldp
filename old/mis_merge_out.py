#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
merge to 
@author: jiedeng
"""
from shared_functions import load_paths, merge, merge_sel, remove_recal_traj_files
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--mode","-m",default = 'sel',type=str,help="input path file")
parser.add_argument("--remove","-r",default = 'n',type=str,help="if input, remove subfolders")


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
    
for path in paths:
    print(path)
    if (not args.mode) and not(args.mode == 'sel'):
        merge(path)  
    else:
        merge_sel(path)

if args.remove:
    print("remove all subfolders")
    for path in paths:
        remove_recal_traj_files(path)
            
    
            

