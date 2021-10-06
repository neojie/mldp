#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 18:48:47 2020
https://stackoverflow.com/questions/52403065/argparse-optional-boolean
@author: jiedeng
"""

from shared_functions import load_paths
import os
import argparse
cwd    = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file.\
                    1-None, current folder\
                    2-all, all subfolders \
                    3-file, folders stored at file")

parser.add_argument("--OUTCAR","-o",type = str, default = 'OUTCAR', help="OUTCAR name")
#parser.add_argument("--verbose","-v",type = bool, default = True, help="OUTCAR name")

parser.add_argument('--verbose',"-v", default=True, action='store_false')
parser.add_argument('--remove_outcar',"-ro", default=False, action='store_true')

args   = parser.parse_args()

if not args.inputpath:
    print("No folders point are provided. Use current working path")
    inputpath = os.path.join(cwd)
    paths = [cwd]
else:
    if args.inputpath == 'all':
        paths = []
        print("check all sub folders at {0}".format(cwd))
        for directory in os.listdir(cwd):
            if os.path.isdir(directory):
                paths.append(directory)
    else:
        print("Check files in {0}  ".format(args.inputpath))
        inputpath = args.inputpath
        paths = load_paths(inputpath)        


def check(path):
    outcar = os.path.join(path,args.OUTCAR)
    fp=open(outcar)
    flag_nelm, flag_nbands, flag_iter, flag_occ = False, False, False, False
    check_iter, check_nbands = False, False
    good_iter, good_nbands = True, False
    for i in range(100000):
        line = fp.readline()
        if 'NELM' in line:#1st appear
            nelm = int(line.split()[2].replace(';','').replace(',','').replace('=',''))
            flag_nelm = True
        if 'number of bands    NBANDS' in line:#2nd appear
            nbands = int(line.split()[-1].replace(';','').replace(',','').replace('=',''))
            flag_nbands = True
        if 'Iteration' in line:#3rd appear
            iteration = int(line.split('(')[-1].split(')')[0])
            flag_iter = True
        if 'band No.' in line:
            for _ in range(nbands):
                line_kpoints = fp.readline()
#            print(line_kpoints)
            occupancy=float(line_kpoints.split()[-1])
            flag_occ = True
        if flag_nelm and flag_iter: 
            check_iter = True
            if nelm == iteration: ## iteration varies from [1 to nelm]
                good_iter = False
                
        if flag_nbands and flag_occ:
            check_nbands = True
            if occupancy<1e-6:
                good_nbands = True
                
        if check_iter and check_nbands:
            if args.verbose:
                print("details")
                print(nelm, iteration)
                print(nbands, occupancy)
            return good_iter, good_nbands
    return good_iter, good_nbands  

for path in paths:
    if os.path.exists(os.path.join(path,args.OUTCAR)):
        good_iter, good_nbands  = check(path)
        if not good_iter:
            print("bad iter",path)
        if not good_nbands:
            print("bad nbands",path)
        if args.remove_outcar and ((not good_iter) or (not good_nbands)):
            os.remove(os.path.join(path,args.OUTCAR))
            print("** remove outcar successfully")
    else:
        if args.verbose:
            print("NO OUTCAR")
  
            
