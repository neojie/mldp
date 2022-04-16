#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 19:53:44 2020

@author: jiedeng
"""

import argparse
import os
import json
import shutil

parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="json")
args   = parser.parse_args()



with open(args.inputpath) as f:
    tmp = json.load(f)
    try:
        paths = tmp['training']['systems']
    except:
        paths = tmp['training']['training_data']['systems']
                

target_folder = '/u/project/ESS/lstixrud/jd848/dat_sum'

if os.path.exists(target_folder):
    shutil.rmtree(target_folder)
deepmd_files = ['energy.npy', 'box.npy',  'coord.npy' , 'force.npy' ,\
                'fparam.npy' , 'virial.npy','type_map.raw', 'type.raw',\
                'set.000', 'set.001', 'set.002', 'set.003', 'set.004', \
                'set.005', 'set.006', 'set.007', 'set.008', 'set.009', \
                'set.010', 'set.011', 'set.012', 'set.013', 'set.014']

def ig_f(dire, files):
    tmpf = []
    for file in files:
        if not (file in deepmd_files):
            tmpf.append(file)

    if len(tmpf) == 0:
        tmpf = [None]
    else:
        print("-"*100)
        print(dire)
        print("remove",tmpf)
    return tmpf

os.mkdir(target_folder)
for path in paths:
    path = os.path.normpath(path) # this is critical
    dire=path.split('/')
    tmp = ''
    for i in dire[1:-1]:
        tmp = tmp +'.'+i
    newfolder = os.path.join(target_folder,tmp[1:])
    os.mkdir(newfolder)
    shutil.copytree(path,os.path.join(newfolder,dire[-1]), ignore = ig_f)
