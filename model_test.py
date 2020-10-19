#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 07:50:22 2020
cd /gpfs/loomis/project/kklee/jd848/pv+hf/5k/60g/r3-80atoms/recal
python ~/script/dpkit/model_test.py -m /gpfs/loomis/project/kklee/jd848/pv+hf/dp-train/345k-pv2/pv.pb

usage: dp test [-h] [-m MODEL] [-s SYSTEM] [-S SET_PREFIX] [-n NUMB_TEST]
               [-r RAND_SEED] [--shuffle-test] [-d DETAIL_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL, --model MODEL
                        Frozen model file to import
  -s SYSTEM, --system SYSTEM
                        The system dir. Recursively detect systems in this
                        directory
  -S SET_PREFIX, --set-prefix SET_PREFIX
                        The set prefix
  -n NUMB_TEST, --numb-test NUMB_TEST
                        The number of data for test
  -r RAND_SEED, --rand-seed RAND_SEED
                        The random seed
  --shuffle-test        Shuffle test data
  -d DETAIL_FILE, --detail-file DETAIL_FILE
                        The file containing details of energy force and virial
                        accuracy

@author: jiedeng
"""

from dp_test import test_ener
from shared_functions import load_paths
import os
import argparse

parser_tst = argparse.ArgumentParser()
parser_tst.add_argument("--inputpath","-ip",help="input path file")

parser_tst.add_argument("-m", "--model", default="frozen_model.pb", type=str, 
                        help="Frozen model file to import")

parser_tst.add_argument("-S", "--set-prefix", default="set", type=str, 
                        help="The set prefix")
parser_tst.add_argument("-n", "--numb-test", default=100, type=int, 
                        help="The number of data for test")
parser_tst.add_argument("-r", "--rand-seed", type=int, 
                        help="The random seed")
parser_tst.add_argument("--shuffle-test", action = 'store_true', 
                        help="Shuffle test data")

parser_tst.add_argument("-d", "--detail-file", type=str, 
                        help="The file containing details of energy force and virial accuracy")

args   = parser_tst.parse_args()
cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,'recal') for path in tmp]
#    paths = tmp
else:
    print("No folders point are provided. Use current folder")
    inputpath = os.path.join(cwd)
    paths = [cwd]
    
inputs = {'model':args.model,
          'rand_seed':0,
          'system': None,
          'set_prefix':'set',
          'shuffle_test': args.shuffle_test,
          'numb_test':args.numb_test,
          'detail_file':args.detail_file}

def extract_sigma_vol_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    get_vol   = False
    get_sigma = False
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            get_sigma = True
        if 'volume of cell' in line:
            vol = float(line.split('volume of cell')[1].split(':')[1].split()[0])
            get_vol = True
        if get_vol and get_sigma:
            outcar_file.close()
            return vol, sigma


out = []
for path in paths:
    
    print(path)
    deepmd_path = os.path.join(path,'deepmd')
    outcar = os.path.join(path,'OUTCAR')
    vol, sigma0 = extract_sigma_vol_outcar(outcar)
    inputs['system'] = deepmd_path
    num,simga,natoms, l2e, l2ea, l2f, l2v = test_ener(inputs)
    print(simga, l2e, l2ea, l2f, l2v)

    if abs(sigma0 - simga)>0.01:
        print("!!!!ERROR ON SIGMA, do not match!!!!")
    out.append([num,vol,simga,natoms, l2e, l2ea, l2f, l2v])
    
import numpy as np
out = np.array(out)
np.savetxt('model_test.out', out)
    

