#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020
shared function, select what to drop
1-analysis deepmd only, analysis outcar only
outcar only, if no model provided


only consider 
@author: jiedeng
"""

print("**"*40)
print("*"*15,' '*19, 'ip format',' '*15, "*"*18 )
print("*"*15,' '*10, 'xx, xx/recal/xx/recal/deepmd*',' '*10,"*"*12 )
print("*"*15,' '*10, '  if no ip, in recal folder  ',' '*10,"*"*12 )

print("**"*40)

import argparse
import os

parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file. support txt file, json file, default is cwd")
parser.add_argument("-m", "--model", type=str,help="Frozen model file to import")
parser.add_argument("-de", "--deepmd", type=str,default = 'deepmd',help="deepmd folder")

parser.add_argument("-S", "--set_prefix", default="set", type=str, 
                        help="The set prefix")
parser.add_argument("-n", "--numb_test", default=100000, type=int,
                        help="The number of data for test") # make default as large as possible
parser.add_argument("-r", "--rand_seed", type=int, 
                        help="The random seed")
parser.add_argument("--shuffle_test", action = 'store_true', 
                        help="Shuffle test data")
parser.add_argument("-d", "--detail_file", type=str, 
                        help="The file containing details of energy force and virial accuracy")

parser.add_argument("-o", "--overwrite", type=True,action='store_false', 
                        help="Overwrite existing test results")
args   = parser.parse_args()
cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    if len(inputpath) >4 and inputpath[-4:] == 'json':
        print("input is json file, load json")
        import json
        with open(inputpath) as f:
            tmp = json.load(f)
            paths = tmp['training']['systems']
    else:
        from shared_functions import load_paths
        paths = load_paths(inputpath,level='recal')
#    paths = [os.path.join(path,args.suffix) for path in tmp]
else:
    print("No folders point are provided. Use current working path")
#    inputpath = os.path.join(cwd)
    paths = [cwd]


inputs = {'model':args.model,
          'rand_seed':0,
          'system': None,
          'set_prefix':'set',
          'shuffle_test': args.shuffle_test,
          'numb_test':args.numb_test,
          'detail_file':args.detail_file}

from dp_test import test_ener, train_ener
for path in paths:
    print('--'*40)
    print(path)
    print('--'*40)
    ## check the source data, some only has deepmd, not deepmd_relax2
    ## and vice versa. deepmd_relax2 is for relax process only
    ## to do the statisitics, we only care where there is deepmd
    try:
        assert os.path.exists(os.path.join(path,args.deepmd))  # deepmd path must exist
        deepmd_path = os.path.join(path,args.deepmd)
    except:
        assert os.path.exists(path) # this gives the fixability 
        deepmd_path = path
        
#    deepmd_path = os.path.join(path,args.deepmd)
    if not args.overwrite and os.path.exists(deepmd_path,'log.'+args.detail_file):
        print('dp test file exist, skip')
    else:
        inputs['system'] = deepmd_path     
        num_test, sigma, natoms, l2e_test, l2ea_test, l2f_test, l2v_test = test_ener(inputs)
        num_tr,   sigma, natoms, l2e_tr,   l2ea_tr,   l2f_tr,   l2v_tr   = train_ener(inputs)
    
    
    ### save log file
        logfile = os.path.join(deepmd_path,'log.'+args.detail_file)
        if os.path.exists(logfile):
            print(logfile + "exist, overwrite")
        log = open(logfile,'w')
        log.write('## model path: ## \n')
        if '/' in args.model:
            log.write('# '+args.model+'\n')
        else:
            log.write('# '+os.path.join(cwd, args.model)+'\n')
        log.writelines(['## system info ## \n','# natoms = {0} \n'.format(natoms), '# sigma = {0} \n'.format(sigma)])
        log.writelines([str(i)+'\n' for i in [natoms,sigma]])        
        log.write ('## train results stored at {0}.*.tr.out##'.format(args.detail_file));log.write ('\n')
        log.write ("# number of train data : %d " % num_tr);log.write ('\n')
        log.write ("# Energy L2err        : %e eV" % l2e_tr);log.write ('\n')
        log.write ("# Energy L2err/Natoms : %e eV" % l2ea_tr);log.write ('\n')
        log.write ("# Force  L2err        : %e eV/A" % l2f_tr);log.write ('\n')
        log.write ("# Virial L2err        : %e eV" % l2v_tr);log.write ('\n')
        log.write ("# Virial L2err/Natoms : %e eV" % (l2v_tr/natoms));log.write ('\n')
        log.writelines([str(i)+'\n' for i in [num_tr,l2e_tr,l2ea_tr,l2f_tr,l2v_tr,l2v_tr/natoms]])
        log.write('## test results stored at {0}.*.out##'.format(args.detail_file));log.write ('\n')
        log.write ("# number of test data : %d " % num_test);log.write ('\n')
        log.write ("# Energy L2err        : %e eV" % l2e_test);log.write ('\n')
        log.write ("# Energy L2err/Natoms : %e eV" % l2ea_test);log.write ('\n')
        log.write ("# Force  L2err        : %e eV/A" % l2f_test);log.write ('\n')
        log.write ("# Virial L2err        : %e eV" % l2v_test);log.write ('\n')
        log.write ("# Virial L2err/Natoms : %e eV" % (l2v_test/natoms));log.write ('\n')
        log.writelines([str(i)+'\n' for i in [num_test,l2e_test,l2ea_test,l2f_test,l2v_test,l2v_test/natoms]])                
        log.close()           
