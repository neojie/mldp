#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
merge to 
@author: jiedeng
"""
from shared_functions import load_paths
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--deepmd","-d",type = str, default = 'deepmd_relax', help="deepmd folder for relax runs")

args   = parser.parse_args()

cwd    = os.getcwd()

deepmds = args.deepmd.split('-')
out = []

print("Check files in {0}  ".format(args.inputpath))
inputpath = args.inputpath
tmp = load_paths(inputpath,level = 'vasp')# actually we can just use `recal` level
if isinstance(tmp[0],list): 
    paths = []
    for i in range(len(tmp[2])):
        if tmp[1][i] >0: # only keep relax > 0
            paths.append(os.path.join(tmp[2][i],'recal'))
else:
    paths = [os.path.join(path,'recal') for path in tmp]

for path in paths:    
    for deepmd in deepmds:
        deepmd_path = os.path.join(path,deepmd)
        if os.path.exists(deepmd_path):
            out.append(deepmd_path)
                
with open('folders_to_train','w') as ft:
    ft.write('        "systems":      '+'[')
    for path in out:
        ft.write('"'+path+'",')
        ft.write("\n")
    ft.write('],')
        
with open('folders_to_train2','w') as ft:
    for path in out:
        ft.write(path)
        ft.write("\n")

