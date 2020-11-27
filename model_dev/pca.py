#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 20:36:32 2020
Given OUTCAR, perform PCA on 'selected' trajectories if id/vid given, if not, on whole trajectory
Given deepmd, perform PCA on 'selected' trajectories if id/vid given, if not, on whole trajectory

1-Allow setting different cutoff coefficients
2-Allow seeting cutoff variance

@author: jiedeng
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import dpdata
import os

parser = argparse.ArgumentParser()
parser.add_argument("--outcar","-o",help="outcar")
parser.add_argument("--deepmd","-d",help="deepmd folder")
parser.add_argument("--coef_threshold","-ct",default=0.01,type=float,help="pca_component_coef_threshold")
parser.add_argument("--variance_threshold","-vt",default=0.99,type=float,help="variance_ratio_threshold")
parser.add_argument("--out_deepmd","-od",default='deepmd_pca',help="output deepmd folder")
parser.add_argument("--vaspidx","-vid",help="vasp id file")
parser.add_argument("--idx","-id",help="id file")
parser.add_argument("--out_idx","-oi",default='pca',help="output idx file prefix")

parser.add_argument('--plot',"-p", default=True, action='store_false',help="plot the results?")
parser.add_argument('--test',"-t", default=True, action='store_false',help="save test as set.001?")

args        = parser.parse_args()

scaler = StandardScaler()
def cal_pca(dat,natoms=160):
    """
    dat in n*160*3 format, n is # of trajectory
    """
    dat = dat.reshape((len(dat),480))
    scaler.fit(dat)
    dat = scaler.transform(dat)    
    pca = PCA().fit(dat.T)
    return pca
    
def sel(pca,variance_ratio_threshold=0.99,pca_component_coef_threshold=0.01):
    variance_ratio = np.cumsum(pca.explained_variance_ratio_)
    min_pca_component = np.argmin(abs(variance_ratio-variance_ratio_threshold))
    pca_component_abs_coef = abs(pca.components_[min_pca_component,:])
    idx1 = np.where(pca_component_abs_coef>pca_component_coef_threshold)[0]
    idx2 = np.where(pca_component_abs_coef<=pca_component_coef_threshold)[0]
    print("In total {0} out of {1} if pca component coef > {2} and variance ratio > {3}".format(
            len(idx1),len(idx1)+len(idx2),pca_component_coef_threshold,variance_ratio_threshold))
    return idx1, idx2, pca_component_abs_coef

def extract_sigma_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            outcar_file.close()
            return sigma
        
if args.outcar:
    ls = dpdata.LabeledSystem(args.outcar,fmt='outcar')
    sigma = extract_sigma_outcar(args.outcar)
    
if args.deepmd:
    print("**** If multiple sets exists, configs are disorded ****")
    ls = dpdata.System(args.deepmd,fmt='deepmd/npy')
    sigma = np.load(os.path.join(args.deepmd,'set.000/fparam.npy'))[0]

idx = np.array(range(len(ls)))
if args.idx:
    print("index file provided")
    idx = np.loadtxt(args.idx).astype('int')
    ls = ls.sub_system(idx)
if args.vaspidx:
    print("vasp index file provided")
    vaspidx=np.loadtxt(args.vaspidx).astype('int')  
    idx = vaspidx - 1        
    ls = ls.sub_system(idx)    

try:
    natoms = ls.get_natoms()
except:
    natoms = 160
    print('assume natoms = 160')
    
pca_ls = cal_pca(ls.data['coords'],natoms);

idx1,idx2,pca_component_abs_coef=sel(pca_ls,
              variance_ratio_threshold = args.variance_threshold,
              pca_component_coef_threshold = args.coef_threshold) 

ls_idx1  = ls.sub_system(idx1)
ls_idx2  = ls.sub_system(idx2)

ls_idx1.to_deepmd_npy(args.out_deepmd,set_size=100000)
if args.test:
    ls_idx2.to_deepmd_npy('pca_test_tmp',set_size=100000)
    import shutil, os
    shutil.copytree('pca_test_tmp/set.000',os.path.join(args.out_deepmd,'set.001'))
    shutil.rmtree('pca_test_tmp')
    
import glob
def build_fparam( folder,sigma): # deepmd/..
    sets=glob.glob(folder+"/set*")
    for seti in sets:
        energy=np.load(seti+'/energy.npy')
        size = energy.size
        all_te = np.ones(size)*sigma
        np.save( os.path.join(seti,'fparam.npy'), all_te)

if args.plot:
    fig,ax = plt.subplots(1,2,figsize=(8,4))
    ax[0].plot(pca_component_abs_coef,'.')
    variance_ratio = np.cumsum(pca_ls.explained_variance_ratio_)
    ax[1].plot(variance_ratio)
    ax[0].set_xlabel('trajectory#')
    ax[0].set_ylabel('coef')
    ax[1].set_xlabel('number of components')
    ax[1].set_ylabel('cumulative explained variance')
    ax[1].grid(True)
    fig.savefig('pca.png')

build_fparam(args.out_deepmd,sigma)    
np.savetxt(args.out_idx+'_pca_id',idx[idx1],fmt='%d')
np.savetxt(args.out_idx+'_pca_vid',idx[idx1]+1,fmt='%d')
