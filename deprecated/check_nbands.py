#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 18:34:32 2020

@author: jiedeng
"""
from shared_functions import load_paths
import os
import glob 
import argparse
cwd    = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath)
else:
    print("No folders point are provided. Use current working path")
    inputpath = os.path.join(cwd)
    paths = [cwd]


def check(fp):
    for i in range(1000):
        line = fp.readline()
        if 'TOO FEW BANDS' in line:
            fp.close()
            return False
    fp.close()
    return True

for path in paths:
    os.chdir(path)
    slurm=glob.glob('slurm*')
    fp=open(slurm[0])
    if check(fp):
        print("ok",path)
    else:
        print("xx",path)
    os.chdir(cwd)


