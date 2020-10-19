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
parser.add_argument("--deepmd","-d",type = str, default = 'deepmd', help="deepmd folder for relax runs")

args   = parser.parse_args()

cwd    = os.getcwd()
out = []
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    if len(tmp)==1:
        out = [os.path.join(os.path.join(path,'recal'),'deepmd') for path in tmp]
    else:
        for i in range(len(tmp)):
#            if tmp[1][i] >0: # only keep relax > 0
#                deepmd = os.path.join(os.path.join(tmp[0][i],'recal'), 'deepmd')
#                deepmd_relax = os.path.join(os.path.join(tmp[-1][i],'recal'), args.deepmd)
#                if os.path.exists(deepmd) and os.path.exists(deepmd_relax):
#                    if deepmd in out:
#                        print("Redundant count!!",deepmd)
#                    else:
#                        out.append(deepmd)
#                    if deepmd_relax in out:
#                        print("Redundant count!!",deepmd_relax)
#                    else:
#                        out.append(deepmd_relax)
#                else:
#                    print(tmp[2][i], " abnormal")
#            else:
            deepmd = os.path.join(os.path.join(tmp[i],'recal'), 'deepmd')
            if os.path.exists(deepmd):
                if deepmd in out:
                    print("Redundant count!!",deepmd)
                else:
                    out.append(deepmd)
            else:
                print(tmp[i], " abnormal")
else:
    out = [cwd]                    
for path in out:
    print(path)
    files = os.listdir(path)
    for file in files:
        if 'xset' in file:
            string = file.replace('xset','x')
            os.rename(os.path.join(path,file),os.path.join(path,string))
            

