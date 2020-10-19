#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 19:14:51 2020

@author: jiedeng
"""
import os
import shutil
cwd    = os.getcwd()
string = 'free  energy   TOTEN'
def check_outcar_done_slow(path):
    """
    slowly check if outcar is done
    require large memory if OUTCAR is large
    """
    outcar = os.path.join(subfolder,'OUTCAR')
    count = 0
    fp  =  open(outcar)
    lines =fp.readlines() 
    for line in lines:
        if string in line:
            count +=1
#            print(count)
            if count>1:
                print(outcar)
                shutil.rmtree(path)
    fp.close()
#    return False

subfolders = [f.path for f in os.scandir(cwd) if f.is_dir()]

for subfolder in subfolders:
    check_outcar_done_slow(subfolder)
