#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:02:21 2019
version 2 
#readlines behave very funny, cannot understand...
# first build folder with XDATCAR, INCAR,POTCAR, bs, KPOINTS


@author: jiedeng
"""

import os


def load_paths(inputpath):
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    fp.close()
    for i in range(len(folders_org)):
        if '#' in folders_org[i] or folders_org[i] == '\n':
            pass
        
        else:
            tmp = folders_org[i].replace('\n','')
            paths.append(tmp.replace('"','').replace(',','').replace('deepmd',''))
    return paths



    
    
#paths = read_path()
#paths = ['/Users/jiedeng/Documents/tmp/jd848/Fe98Si10O20/melt4/r4-4000-FeSiC/']
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd,'folders_pv')
    
#paths = ['/gpfs/loomis/project/kklee/jd848/pv+hf/4k/100g/r3']
### 

import shutil
import glob
paths = load_paths(inputpath)
for path in paths:
    # get basic info

    if os.path.exists(os.path.join(path,'recal')):
        target_path = os.path.join(path,'recal')
        
        files =  glob.glob(os.path.join(target_path,'*.tsv'))
        
        if len(files) ==0:
            shutil.rmtree(os.path.join(path,'recal'))
            print("remove",target_path)
#        path = os.path.join(path,'recal')
#        pass 
    else:
        pass


