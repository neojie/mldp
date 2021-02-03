#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 09:43:08 2020

analyze post stat_model.py run

1- stat_model.py generate test and train results in each folder
2- summarize all the data, and store it as train test set
3- analyze the results, make a prediction vs. results
4-

TODO : consider vol and natoms change frame to frame
@author: jiedeng
"""

import argparse


parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file. support txt file, json file, default is cwd")
parser.add_argument("-de", "--deepmd", type=str,default = 'deepmd',help="deepmd folder")
parser.add_argument("-d", "--detail_file", type=str, 
                        help="The file containing details of energy force and virial accuracy")

parser.add_argument("-mg", "--merge", default=False,action = 'store_true',
                        help="merge test and set, if test and set is the same, only merge once > should be considered in stat_model.py, not here")
args   = parser.parse_args()

import numpy as np
import os
import matplotlib.pyplot as plt
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
    
eV_A3_2_GPa  = 160.21766208 # 1 eV/Å3 = 160.2176621 GPa

count = 0
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
    logfile  = os.path.join(deepmd_path,'log.'+args.detail_file)
    # deepmd must at least have set.000 and type.raw
    natoms = np.loadtxt(os.path.join(deepmd_path,'type.raw')).shape[0]  #assume natoms does not change
    boxs = np.load(os.path.join(deepmd_path,'set.000/box.npy')) #box may shape may vary from frame to frame => NPT, MPMT
    box = boxs[0]
    a = box[:3]
    b = box[3:6]
    c = box[6:]
    vol = np.dot(c,np.cross(a,b))

    e_tr_file  =  os.path.join(deepmd_path,args.detail_file+".e.tr.out")
    f_tr_file  =  os.path.join(deepmd_path,args.detail_file+".f.tr.out")
    v_tr_file  =  os.path.join(deepmd_path,args.detail_file+".v.tr.out")
    
    e_test_file = os.path.join(deepmd_path, args.detail_file+".e.out")
    f_test_file = os.path.join(deepmd_path, args.detail_file+".f.out")
    v_test_file = os.path.join(deepmd_path, args.detail_file+".v.out")
    
    train_exist = False
    if os.path.exists(e_tr_file):
        train_exist = True
    else:
        print("WARNING: train results not availabel")
    test_exist = False
    if os.path.exists(e_test_file):
        test_exist = True   
    else:
        print("WARNING: test results not availabel")
        
    if train_exist:
        e_tr = np.loadtxt(e_tr_file)/natoms
        f_tr = np.loadtxt(f_tr_file)
        v_tr = np.loadtxt(v_tr_file);     
        if e_tr.ndim == 1:
            e_tr = np.reshape(e_tr,(1,len(e_tr)))
        if v_tr.ndim == 1:
            v_tr = np.reshape(v_tr,(1,len(v_tr)))
        v_gpa_tr = v_tr/vol*eV_A3_2_GPa 


    if test_exist:        
        e_test = np.loadtxt(e_test_file)/natoms
        f_test = np.loadtxt(f_test_file)
        v_test = np.loadtxt(v_test_file); 
        
        if e_test.ndim == 1:
            e_test = np.reshape(e_test,(1,len(e_test)))
        if v_test.ndim == 1:
            v_test = np.reshape(v_test,(1,len(v_test)))
        # force cannot be 1 d datal as it is for every atom
        v_gpa_test = v_test/vol*eV_A3_2_GPa            
 
    
#    log = np.loadtxt(logfile)
    
    if count == 0:
        if train_exist:
            e_tr_all = e_tr
            f_tr_all = f_tr
            v_tr_all = v_tr; v_gpa_tr_all = v_gpa_tr
        if test_exist:
            e_test_all = e_test
            f_test_all = f_test
            v_test_all = v_test; v_gpa_test_all = v_gpa_test
    else:
        if train_exist: #  corner case not considered: 1st directory only has train or test
            e_tr_all = np.concatenate((e_tr_all,e_tr),axis=0)
            f_tr_all = np.concatenate((f_tr_all,f_tr),axis=0)
            v_tr_all = np.concatenate((v_tr_all,v_tr),axis=0); v_gpa_tr_all = np.concatenate((v_gpa_tr_all,v_gpa_tr),axis=0)           
        if test_exist:
            e_test_all = np.concatenate((e_test_all,e_test),axis=0)
            f_test_all = np.concatenate((f_test_all,f_test),axis=0)
            v_test_all = np.concatenate((v_test_all,v_test),axis=0); v_gpa_test_all = np.concatenate((v_gpa_test_all,v_gpa_test),axis=0)           
    
    count += 1

def rmse(org,pred): # same as dp_test l2err
    dif = pred - org
    return np.sqrt(np.mean(dif**2))

def min_max(dat):
    MIN = np.min(np.min(dat,axis=0))
    MAX = np.max(np.max(dat,axis=0))
    return MIN,MAX

def save(dat,name,header):
    np.savetxt(os.path.join(cwd,args.detail_file+name), dat, header = header)

def plot(e,f,v):
    fig, ax = plt.subplots(1,3,figsize=(14,4))
    ax[0].plot(e[:,0],e[:,1],'.')
    ax[0].plot(min_max(e), min_max(e),'k--')
    ax[0].set_xlabel('DFT energy/atom (eV)')
    ax[0].set_ylabel('NN energy/atom (eV)')
    
    
    ax[1].plot(f[:,0],f[:,3],'.')
    ax[1].plot(f[:,1],f[:,4],'.')
    ax[1].plot(f[:,2],f[:,5],'.')
    ax[1].plot(min_max(f_tr_all), min_max(f_tr_all),'k--')
    ax[1].set_xlabel('DFT force (eV/A)')
    ax[1].set_ylabel('NN force (eV/A)')
    
    ax[2].plot(v[:,0],v[:,9],'.')
    ax[2].plot(v[:,1],v[:,10],'.')
    ax[2].plot(v[:,2],v[:,11],'.')
    ax[2].plot(v[:,3],v[:,12],'.')
    ax[2].plot(v[:,4],v[:,13],'.')
    ax[2].plot(v[:,5],v[:,14],'.')
    ax[2].plot(v[:,6],v[:,15],'.')
    ax[2].plot(v[:,7],v[:,16],'.')
    ax[2].plot(v[:,8],v[:,17],'.')
    ax[2].plot(min_max(v), min_max(v),'k--')
    ax[2].set_xlabel('DFT stress (GPa)')
    ax[2].set_ylabel('NN stress (GPa)')
    plt.show()

if args.merge:
    print("merge train and test")
    e_all = np.concatenate((e_test_all,e_test_all),axis=0)
    f_all = np.concatenate((f_test_all,f_tr_all),axis=0)
    v_all = np.concatenate((v_test_all,v_tr_all),axis=0)
    v_gpa_all = np.concatenate((v_gpa_test_all,v_gpa_tr_all),axis=0)
    
    print("N frames: ", len(e_all))
    print("e/atom rmse:", rmse(e_all[:,0],e_all[:,1]))
    print("f rmse:", rmse(f_all[:,:3],f_all[:,3:]))
    print("v rmse:", rmse(v_all[:,:9],v_all[:,9:]))
    print("v (gpa) rmse:", rmse(v_gpa_all[:,:9],v_gpa_all[:,9:]))

else:   
    print("tr frames: ", len(e_tr_all))
    print("tr e/atom rmse:", rmse(e_tr_all[:,0],e_tr_all[:,1]))
    print("tr f rmse:", rmse(f_tr_all[:,:3],f_tr_all[:,3:]))
    print("tr v rmse:", rmse(v_tr_all[:,:9],v_tr_all[:,9:]))
    print("tr v (gpa) rmse:", rmse(v_gpa_tr_all[:,:9],v_gpa_tr_all[:,9:]))
    print("test frames: ", len(e_test_all))    
    print("test e/atom rmse:", rmse(e_test_all[:,0],e_test_all[:,1]))
    print("test f rmse:", rmse(f_test_all[:,:3],f_test_all[:,3:]))
    print("test v rmse:", rmse(v_test_all[:,:9],v_test_all[:,9:]))
    print("test v (gpa) rmse:", rmse(v_gpa_test_all[:,:9],v_gpa_test_all[:,9:]))

# save 
print("##"*100)
if args.merge:
    save(e_all,".merge.e.out",'data_e pred_e')
    save(f_all,".merge.f.out",'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
    save(v_all,".merge.v.out",'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')
else:
    save(e_tr_all,".all.e.tr.out",'data_e pred_e')
    save(f_tr_all,".all.f.tr.out",'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
    save(v_tr_all,".all.v.tr.out",'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')
    
    save(e_test_all,".all.e.test.out",'data_e pred_e')
    save(f_test_all,".all.f.test.out",'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
    save(v_test_all,".all.v.test.out",'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')

if args.merge:
    plot(e_all,f_all,v_gpa_all)
else:
    plot(e_tr_all,f_tr_all,v_gpa_tr_all)
    plot(e_test_all,f_test_all,v_gpa_test_all) 
