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

args   = parser.parse_args()

import numpy as np
import os
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

    e_tr = np.loadtxt(e_tr_file)/natoms
    f_tr = np.loadtxt(f_tr_file)
    v_tr = np.loadtxt(v_tr_file);     v_gpa_tr = v_tr/vol*eV_A3_2_GPa 

    e_test = np.loadtxt(e_test_file)/natoms
    f_test = np.loadtxt(f_test_file)
    v_test = np.loadtxt(v_test_file); v_gpa_test = v_test/vol*eV_A3_2_GPa 
    
    if e_tr.ndim == 1:
        e_tr = np.reshape(e_tr,(1,len(e_tr)))
    if e_test.ndim == 1:
        e_test = np.reshape(e_test,(1,len(e_test)))

    if f_tr.ndim == 1:
        e_tr = np.reshape(f_tr,(1,len(f_tr)))
    if f_test.ndim == 1:
        e_test = np.reshape(f_test,(1,len(f_test)))
    # force cannot be 1 d datal as it is for every atom
    log = np.loadtxt(logfile)
    
    if count == 0:
        e_tr_all = e_tr
        f_tr_all = f_tr
        v_tr_all = v_tr; v_gpa_tr_all = v_gpa_tr

        e_test_all = e_test
        f_test_all = f_test
        v_test_all = v_test; v_gpa_test_all = v_gpa_test
    else:
        e_tr_all = np.concatenate((e_tr_all,e_tr),axis=0)
        f_tr_all = np.concatenate((f_tr_all,f_tr),axis=0)
        v_tr_all = np.concatenate((v_tr_all,v_tr),axis=0); v_gpa_tr_all = np.concatenate((v_gpa_tr_all,v_gpa_tr),axis=0)           

        print(e_test_all)
        print(e_test)
        e_test_all = np.concatenate((e_test_all,e_test),axis=0)
        f_test_all = np.concatenate((f_test_all,f_test),axis=0)
        v_test_all = np.concatenate((v_test_all,v_test),axis=0); v_gpa_test_all = np.concatenate((v_gpa_test_all,v_gpa_test),axis=0)           
    
    count += 1

def rmse(org,pred): # same as dp_test l2err
    dif = pred - org
    return np.sqrt(np.mean(dif**2))

    
print("tr e/atom rmse:", rmse(e_tr_all[:,0],e_tr_all[:,1]))
print("tr f rmse:", rmse(f_tr_all[:,:3],f_tr_all[:,3:]))
print("tr v rmse:", rmse(v_tr_all[:,:9],v_tr_all[:,9:]))
print("tr v (gpa) rmse:", rmse(v_gpa_tr_all[:,:9],v_gpa_tr_all[:,9:]))

print("test e/atom rmse:", rmse(e_test_all[:,0],e_test_all[:,1]))
print("test f rmse:", rmse(f_test_all[:,:3],f_test_all[:,3:]))
print("test v rmse:", rmse(v_test_all[:,:9],v_test_all[:,9:]))
print("test v (gpa) rmse:", rmse(v_gpa_test_all[:,:9],v_gpa_test_all[:,9:]))

# save 
print("##"*100)
print("save merged e f v results")

np.savetxt(os.path.join(cwd,args.detail_file+".all.e.tr.out"), e_tr_all, 
           header = 'data_e pred_e')
np.savetxt(os.path.join(cwd,args.detail_file+".all.f.tr.out"), f_tr_all,
           header = 'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
np.savetxt(os.path.join(cwd,args.detail_file+".all.v.tr.out"), v_tr_all,
           header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        
np.savetxt(os.path.join(cwd,args.detail_file+".all_gpa.v.tr.out"), v_gpa_tr_all,
           header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        

np.savetxt(os.path.join(cwd,args.detail_file+".all.e.test.out"), e_test_all, 
           header = 'data_e pred_e')
np.savetxt(os.path.join(cwd,args.detail_file+".all.f.test.out"), f_test_all,
           header = 'data_fx data_fy data_fz pred_fx pred_fy pred_fz')
np.savetxt(os.path.join(cwd,args.detail_file+".all.v.test.out"), v_test_all,
           header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        
np.savetxt(os.path.join(cwd,args.detail_file+".all_gpa.v.test.out"), v_gpa_test_all,
           header = 'data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz')        

def min_max(dat):
    MIN = np.min(np.min(dat,axis=0))
    MAX = np.max(np.max(dat,axis=0))
    return MIN,MAX

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].plot(e_tr_all[:,0],e_tr_all[:,1],'.')
ax[0].plot(min_max(e_tr_all), min_max(e_tr_all),'k--')
ax[0].set_xlabel('DFT energy/atom (eV)')
ax[0].set_ylabel('NN energy/atom (eV)')


ax[1].plot(f_tr_all[:,0],f_tr_all[:,3],'.')
ax[1].plot(f_tr_all[:,1],f_tr_all[:,4],'.')
ax[1].plot(f_tr_all[:,2],f_tr_all[:,5],'.')
ax[1].plot(min_max(f_tr_all), min_max(f_tr_all),'k--')
ax[1].set_xlabel('DFT force (eV/A)')
ax[1].set_ylabel('NN force (eV/A)')

ax[2].plot(v_gpa_tr_all[:,0],v_gpa_tr_all[:,9],'.')
ax[2].plot(v_gpa_tr_all[:,1],v_gpa_tr_all[:,10],'.')
ax[2].plot(v_gpa_tr_all[:,2],v_gpa_tr_all[:,11],'.')
ax[2].plot(v_gpa_tr_all[:,3],v_gpa_tr_all[:,12],'.')
ax[2].plot(v_gpa_tr_all[:,4],v_gpa_tr_all[:,13],'.')
ax[2].plot(v_gpa_tr_all[:,5],v_gpa_tr_all[:,14],'.')
ax[2].plot(v_gpa_tr_all[:,6],v_gpa_tr_all[:,15],'.')
ax[2].plot(v_gpa_tr_all[:,7],v_gpa_tr_all[:,16],'.')
ax[2].plot(v_gpa_tr_all[:,8],v_gpa_tr_all[:,17],'.')
ax[2].plot(min_max(v_gpa_tr_all), min_max(v_gpa_tr_all),'k--')
ax[2].set_xlabel('DFT stress (GPa)')
ax[2].set_ylabel('NN stress (GPa)')
plt.show()

fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].plot(e_test_all[:,0],e_test_all[:,1],'.')
ax[0].plot(min_max(e_test_all), min_max(e_test_all),'k--')
ax[0].set_xlabel('DFT energy/atom (eV)')
ax[0].set_ylabel('NN energy/atom (eV)')

ax[1].plot(f_test_all[:,0],f_test_all[:,3],'.')
ax[1].plot(f_test_all[:,1],f_test_all[:,4],'.')
ax[1].plot(f_test_all[:,2],f_test_all[:,5],'.')
ax[1].plot(min_max(f_test_all), min_max(f_test_all),'k--')
ax[1].set_xlabel('DFT force (eV/A)')
ax[1].set_ylabel('NN force (eV/A)')

ax[2].plot(v_gpa_test_all[:,0],v_gpa_test_all[:,9],'.')
ax[2].plot(v_gpa_test_all[:,1],v_gpa_test_all[:,10],'.')
ax[2].plot(v_gpa_test_all[:,2],v_gpa_test_all[:,11],'.')
ax[2].plot(v_gpa_test_all[:,3],v_gpa_test_all[:,12],'.')
ax[2].plot(v_gpa_test_all[:,4],v_gpa_test_all[:,13],'.')
ax[2].plot(v_gpa_test_all[:,5],v_gpa_test_all[:,14],'.')
ax[2].plot(v_gpa_test_all[:,6],v_gpa_test_all[:,15],'.')
ax[2].plot(v_gpa_test_all[:,7],v_gpa_test_all[:,16],'.')
ax[2].plot(v_gpa_test_all[:,8],v_gpa_test_all[:,17],'.')
ax[2].plot(min_max(v_gpa_test_all), min_max(v_gpa_test_all),'k--')
ax[2].set_xlabel('DFT stress (GPa)')
ax[2].set_ylabel('NN stress (GPa)')
plt.show()
