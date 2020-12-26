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
#import numpy as np
import os
from shared_functions import load_paths


parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("-m", "--model", type=str,help="Frozen model file to import")
#parser.add_argument("-mo", "--mo", type=str,help="model only")
#parser.add_argument("-de", "--deepmd", type=str,default = 'deepmd',help="deepmd folder, could be deepmd-deepmd_relax, use - as separator")
parser.add_argument("-de", "--deepmd", type=str,default = 'deepmd',help="deepmd folder")


#parser.add_argument("--outcar","-o",type=str,default = 'OUTCAR',help="name of outcar file")
#parser.add_argument("--excudedfraction","-e",type=float,default = 0.2,help="excluded fraction, first few may not reach equil")
#parser.add_argument("--suffix","-s",type=str,default = '',help="add recal or other suffix in the path")

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

args   = parser.parse_args()
cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath,level='recal')
#    paths = [os.path.join(path,args.suffix) for path in tmp]
else:
    print("No folders point are provided. Use current working path")
    inputpath = os.path.join(cwd)
    paths = [cwd]

#out = []
#phases = {'liquid':0,'pv':1,'ppv':2,'akimotoite':3,'enstatite':4, 'major':5, 'unkown':9}

#def get_phase(path):
#    if 'ppv' in path:
#        phase = 'ppv'
#    elif 'pv' in path:
#        phase = 'pv'
#    elif 'enstatite' in path:
#        phase = 'enstatite'
#    elif 'major' in path:
#        phase = 'major'
#    elif 'akimotoite' in path:
#        phase = 'akimotoite'   
#    else:
#        print('phase unidentifiable')
#        phase = 'unkown'
#    return phase


if args.model:  ##
    from dp_test import test_ener, train_ener
    inputs = {'model':args.model,
              'rand_seed':0,
              'system': None,
              'set_prefix':'set',
              'shuffle_test': args.shuffle_test,
              'numb_test':args.numb_test,
              'detail_file':args.detail_file}




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
        assert os.path.exists(path)
        deepmd_path = path
        
#    deepmd_path = os.path.join(path,args.deepmd)
    inputs['system'] = deepmd_path        
    num_test, sigma, natoms, l2e_test, l2ea_test, l2f_test, l2v_test = test_ener(inputs)
    num_tr,   sigma, natoms, l2e_tr,   l2ea_tr,   l2f_tr,   l2v_tr   = train_ener(inputs)


### save log file
    logfile= 'log.'+args.detail_file
    assert (not os.path.exists(os.path.join(deepmd_path,logfile)))
    log = open(logfile,'w')
    log.writelines(['## system info ##','# natoms = {0}'.format(natoms), '# sigma = {0}'.format(sigma)])
    log.writelines([str(i) for i in [natoms,sigma]])        
    log.write ('## train results stored at {0}.*.tr.out##'.format(args.detail_file))
    log.write ("# number of test data : %d " % num_tr)
    log.write ("# Energy L2err        : %e eV" % l2e_tr)
    log.write ("# Energy L2err/Natoms : %e eV" % l2ea_tr)
    log.write ("# Force  L2err        : %e eV/A" % l2f_tr)
    log.write ("# Virial L2err        : %e eV" % l2v_tr)
    log.write ("# Virial L2err/Natoms : %e eV" % l2v_tr/natoms)
    log.writelines([str(i) for i in [num_tr,l2e_tr,l2ea_tr,l2f_tr,l2v_tr,l2v_tr/natoms]])
    log.write('## train results stored at {0}.*.out##'.format(args.detail_file))
    log.write ("# number of test data : %d " % num_test)
    log.write ("# Energy L2err        : %e eV" % l2e_test)
    log.write ("# Energy L2err/Natoms : %e eV" % l2ea_test)
    log.write ("# Force  L2err        : %e eV/A" % l2f_test)
    log.write ("# Virial L2err        : %e eV" % l2v_test)
    log.write ("# Virial L2err/Natoms : %e eV" % l2v_test/natoms)  
    log.writelines([str(i) for i in [num_test,l2e_test,l2ea_test,l2f_test,l2v_test,l2v_test/natoms]])                
    log.close()           
