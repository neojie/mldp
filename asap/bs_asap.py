#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 09:15:25 2020

@author: jiedeng
"""

import pandas as pd
#import shutil
import os


def asap(path,stride=1):
    asap_dir =  os.path.join(path,'asap')
    try:
        os.mkdir(asap_dir)
    except:
        pass
    os.chdir(asap_dir)
    call("asap gen_desc -s {0} --fxyz ../OUTCAR soap -e -c 6 -n 6 -l 6 -g 0.44 --crossover".format(stride), shell=True)

def check_outcar_in_recal_vasp(deepmds):
    new = []
    for deepmd in deepmds:
        recal_dir = os.path.abspath(os.path.join(deepmd, os.pardir))
        if os.path.exists(os.path.join(recal_dir,'OUTCAR')):
            pass
        else:
            print("??",recal_dir,"no outcar in recal")
        pardir = os.path.abspath(os.path.join(recal_dir, os.pardir)) # https://stackoverflow.com/questions/2860153/how-do-i-get-the-parent-directory-in-python
        if os.path.exists(os.path.join(pardir,'OUTCAR')):
            pass
        else:
            print("??",pardir,"no outcar in pardir")
        if os.path.exists(os.path.join(recal_dir,'OUTCAR')) and os.path.exists(os.path.join(pardir,'OUTCAR')):
            new.append(deepmd)
    return deepmds
            
#out_pd = pd.read_excel('/Users/jiedeng/GD/ppt/2020/extreme_filter4_with_state.xlsx')
ak = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'ak')
deepmds = ak['local'].values

## check OUTCAR exists

deepmds = check_outcar_in_recal_vasp(deepmds)

stride = 2
override = True
from subprocess import call
for deepmd in deepmds:
    print('--'*30)
    print(deepmd)
    recal_dir = os.path.abspath(os.path.join(deepmd, os.pardir))
    vaspdir = os.path.abspath(os.path.join(recal_dir, os.pardir)) # https://stackoverflow.com/questions/2860153/how-do-i-get-the-parent-directory-in-python

    if override:
        asap(recal_dir) 
        asap(vaspdir,stride =2)
    else:
        if os.path.exists(os.path.join(recal_dir,'asap')):
            print('exist, skip')
        else:
            asap(recal_dir) 
        if os.path.exists(os.path.join(vaspdir,'asap')):
            print('exist, skip')
        else:
            asap(vaspdir,stride =2)        
#    if not os.path.exists(os.path.join(recal_dir,'asap')):
#        asap(recal_dir)      
#    
#    if not os.path.exists(os.path.join(vasodir,'asap')):
#        asap(recal_dir,stride =2)
        
    