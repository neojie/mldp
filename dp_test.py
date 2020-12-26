#!/usr/bin/env python3

#import re
#import os
#import sys
#import argparse
import numpy as np

from deepmd.Data import DataSets
#from deepmd.Data import DeepmdData
#from deepmd import DeepEval
from deepmd import DeepPot
#from deepmd import DeepDipole
#from deepmd import DeepPolar
#from deepmd import DeepWFC
#from tensorflow.python.framework import ops
import os

def l2err (diff) :    
    return np.sqrt(np.average (diff*diff))

def test_ener (args) :
    """
    modify based on args file
    """
    if args['rand_seed'] is not None :
        np.random.seed(args['rand_seed'] % (2**32))

    data = DataSets (args['system'], args['set_prefix'], shuffle_test = args['shuffle_test'])
    test_data = data.get_test ()
    numb_test = args['numb_test']
    natoms    = len(test_data["type"][0])
    nframes   = test_data["box"].shape[0]
    numb_test = min(nframes, numb_test)
    dp        = DeepPot(args['model'])
    coord     = test_data["coord"][:numb_test].reshape([numb_test, -1])
    box       = test_data["box"][:numb_test]
    atype     = test_data["type"][0]
    if dp.get_dim_fparam() > 0:
        fparam = test_data["fparam"][:numb_test] 
    else :
        fparam = None
    if dp.get_dim_aparam() > 0:
        aparam = test_data["aparam"][:numb_test] 
    else :
        aparam = None
    detail_file = args['detail_file']
    if detail_file is not None:
        atomic = True
    else:
        atomic = False

    ret = dp.eval(coord, box, atype, fparam = fparam, aparam = aparam, atomic = atomic)
    energy = ret[0]
    force  = ret[1]
    virial = ret[2]
    energy = energy.reshape([numb_test,1])
    force  = force.reshape([numb_test,-1])
    virial = virial.reshape([numb_test,9])
    if atomic:
        ae = ret[3]
        av = ret[4]
        ae = ae.reshape([numb_test,-1])
        av = av.reshape([numb_test,-1])

    l2e = (l2err (energy - test_data["energy"][:numb_test].reshape([-1,1])))
    l2f = (l2err (force  - test_data["force"] [:numb_test]))
    l2v = (l2err (virial - test_data["virial"][:numb_test]))
    l2ea = l2e/natoms
    l2va = l2v/natoms

    # print ("# energies: %s" % energy)
    print ("# number of test data : %d " % numb_test)
    print ("Energy L2err        : %e eV" % l2e)
    print ("Energy L2err/Natoms : %e eV" % l2ea)
    print ("Force  L2err        : %e eV/A" % l2f)
    print ("Virial L2err        : %e eV" % l2v)
    print ("Virial L2err/Natoms : %e eV" % l2va)

    if detail_file is not None :
        pe = np.concatenate((np.reshape(test_data["energy"][:numb_test], [-1,1]),
                             np.reshape(energy, [-1,1])), 
                            axis = 1)
        np.savetxt(os.path.join(args['system'], detail_file+".e.out"), pe, 
                   header = 'data_e pred_e')
        pf = np.concatenate((np.reshape(test_data["force"] [:numb_test], [-1,3]), 
                             np.reshape(force,  [-1,3])), 
                            axis = 1)
        np.savetxt(os.path.join(args['system'], detail_file+".f.out"), pf,
                   header = 'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
        pv = np.concatenate((np.reshape(test_data["virial"][:numb_test], [-1,9]), 
                             np.reshape(virial, [-1,9])), 
                            axis = 1)
        np.savetxt(os.path.join(args['system'], detail_file+".v.out"), pv,
                   header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        
    return numb_test,fparam[0][0],natoms, l2e, l2ea, l2f, l2v

def get_train_data(data):
    """
    data: DataSets instance
    ------
    DataSets methods can only load each set directory seperately
    This is a wrapper to concatenate all training data together
    """
    train_dirs = data.train_dirs
    data.load_batch_set(train_dirs[0]);  out = data.batch_set
    
    if len(train_dirs)>1:
        for i in range(1,len(train_dirs)):
            data.load_batch_set(train_dirs[i]); add  = data.batch_set
            for key in out:
                out[key] = np.concatenate((out[key],add[key]))
    return out

def train_ener (inputs) :
    """
    deepmd-kit has function test_ener which deal with test_data only
    `train_ener` are for train data only
    """
    
    if inputs['rand_seed'] is not None :
        np.random.seed(inputs['rand_seed'] % (2**32))

    data = DataSets (inputs['system'], inputs['set_prefix'], shuffle_test = inputs['shuffle_test'])

    
    train_data = get_train_data (data)
    

    numb_test = data.get_sys_numb_batch(1)  ## use 1 batch, # of batches are the numb of train
    natoms = len(train_data["type"][0])
    nframes = train_data["box"].shape[0]
    numb_test = min(nframes, numb_test)
    dp = DeepPot(inputs['model'])
    coord = train_data["coord"].reshape([numb_test, -1])
    box   = train_data["box"]
    atype = train_data["type"][0]
    if dp.get_dim_fparam() > 0:
        fparam = train_data["fparam"]
    else :
        fparam = None
    if dp.get_dim_aparam() > 0:
        aparam = train_data["aparam"]
    else :
        aparam = None
    detail_file = inputs['detail_file']
    if detail_file is not None:
        atomic = True
    else:
        atomic = False

    ret = dp.eval(coord, box, atype, fparam = fparam, aparam = aparam, atomic = atomic)
    energy = ret[0]
    force  = ret[1]
    virial = ret[2]
    energy = energy.reshape([numb_test,1])
    force  = force.reshape([numb_test,-1])
    virial = virial.reshape([numb_test,9])
    if atomic:
        ae = ret[3]
        av = ret[4]
        ae = ae.reshape([numb_test,-1])
        av = av.reshape([numb_test,-1])

    l2e = (l2err (energy - train_data["energy"].reshape([-1,1])))
    l2f = (l2err (force  - train_data["force"]))
    l2v = (l2err (virial - train_data["virial"]))
    l2ea= l2e/natoms
    l2va= l2v/natoms

    # print ("# energies: %s" % energy)
    print ("# number of train data : %d " % numb_test)
    print ("Energy L2err        : %e eV" % l2e)
    print ("Energy L2err/Natoms : %e eV" % l2ea)
    print ("Force  L2err        : %e eV/A" % l2f)
    print ("Virial L2err        : %e eV" % l2v)
    print ("Virial L2err/Natoms : %e eV" % l2va)

    if detail_file is not None :
        pe = np.concatenate((np.reshape(train_data["energy"], [-1,1]),
                             np.reshape(energy, [-1,1])), 
                            axis = 1)
        np.savetxt(os.path.join(inputs['system'],detail_file+".e.tr.out"), pe, 
                   header = 'data_e pred_e')
        pf = np.concatenate((np.reshape(train_data["force"], [-1,3]), 
                             np.reshape(force,  [-1,3])), 
                            axis = 1)
        np.savetxt(os.path.join(inputs['system'],detail_file+".f.tr.out"), pf,
                   header = 'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
        pv = np.concatenate((np.reshape(train_data["virial"], [-1,9]), 
                             np.reshape(virial, [-1,9])), 
                            axis = 1)
        np.savetxt(os.path.join(inputs['system'],detail_file+".v.tr.out"), pv,
                   header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        
    return numb_test,fparam[0][0],natoms, l2e, l2ea, l2f, l2v

