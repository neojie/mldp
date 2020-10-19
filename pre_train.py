#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 09:24:03 2020
prepare train folder
@author: jiedeng
"""


from shared_functions import load_paths
import os
import argparse

# check if exist 

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()


cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,'recal') for path in tmp]
else:
    print("No folders  are provided !")
    raise ValueError()

with open('folders_to_train','w') as ft:
    ft.write('        "systems":      '+'[')
    for path in paths:
        if os.path.exists(os.path.join(path,'deepmd')):
            ft.write('"'+os.path.join(path,'deepmd')+'",')
            ft.write("\n")
        else:
            print("{0} provided but there is no deepmd folder".format(path))
    ft.write('],')
