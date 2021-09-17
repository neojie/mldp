#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 10:48:49 2020
analyze test results 
TWO modes available:
    1) compare nn and VASP   : MUST specify recal folder where there is OUTCAR file
    2) model deviation of nn : MUST NOT specify recal folder
For 1), if multiple nn specified, their averaged is used
@author: jiedeng
"""

import numpy as np
import os
from model_dev_funcs import extract_org_nn_pred, dev_nn, dev_vasp_nn,extract_nn_pred
import argparse
import matplotlib.pyplot as plt

description="""
        example: analysis.py -tf recal -rf recal -mp mm4 mm5 
        extract average recal/mm4.e.out recal/mm5.e.out, the 1st col should be vasp values
        the 2nd col should be model predicted values. 
        mm4 and mm5 are averaged and compared against the vasp values.
        idx between given thresholds are extracted
        \n\n
        example: analysis.py -tf recal -mp mm4 mm5 
        extract recal/mm4.e.out recal/mm5.e.out, the 1st col is ignored
        the 2nd col should be model predicted values, mm4 and mm5 are caompared
        model devation is computed.
        idx between force thresholds are extracted
        \n\n
        Note these two modes, thresholds have different meanings
        """
print(description)
parser = argparse.ArgumentParser()

parser.add_argument("--test_folder","-tf",help="folder storing model test ")
parser.add_argument("--recal_foler","-rf",help="recal folder, the model deivation will be calcualted if this is given")
parser.add_argument("--natoms","-n",default=160,type=int,help="# of atoms in the system")
parser.add_argument("--model_prefix","-mp",nargs="+",default='mm2',help="recal folder, support one model or multiple models, e.g., re4 re5 re6")
parser.add_argument("--energy_lower_cutoff","-elc",default=0.0044*2,type=float,help="lower cutoff for energy")
parser.add_argument("--energy_upper_cutoff","-euc",default=0.1,type=float,help="upper cutoff for energy")
parser.add_argument("--force_lower_cutoff","-flc",default=0.27,type=float,help="lower cutoff for force")
parser.add_argument("--force_upper_cutoff","-fuc",default=1,type=float,help="upper cutoff for force")
args        = parser.parse_args()

mode = 'nn_only'
if args.recal_foler:
    recal_foler = args.recal_foler
    mode = 'nn_vasp'
    print("mode is {0}, nn average and vasp are compared".format(mode))
else:
    ## only nn test results are given
    print("mode is {0}, nn model deviation is analyzed".format(mode))
test_folder = args.test_folder

prefixs     = args.model_prefix
natoms      = args.natoms



if mode == 'nn_vasp':
    # nn and vasp results are stored togheter, no special nsw chosen needed
    es_org,fs_org,vs_org, es,fs,vs   = extract_org_nn_pred(test_folder,prefixs,natoms=natoms)
    etot,forces ,stress              = es_org[0],fs_org[0],vs[0]
    vasp                             = [etot,forces,stress]    
    nsw                              = np.array(range(len(es_org[0])))
else:
    es,fs,vs               = extract_nn_pred(test_folder,prefixs,natoms=natoms)
    
nn   = [np.mean(es,axis=0),np.mean(fs,axis=0),np.mean(vs,axis=0)]


if mode == 'nn_vasp':
    e_diff_per_atom      = (etot[nsw] - nn[0])/natoms
    rmsd_e,rmsd_f,rmsd_v = dev_vasp_nn(vasp,nn,natoms=natoms)

    e_idx1 = np.where(abs(e_diff_per_atom)<args.energy_upper_cutoff)[0]
    e_idx2 = np.where(abs(e_diff_per_atom)>args.energy_lower_cutoff)[0]
    f_idx1 = np.where(abs(rmsd_f)<args.force_upper_cutoff)[0]
    f_idx2 = np.where(abs(rmsd_f)>args.force_lower_cutoff)[0]

    e_idx=np.intersect1d(e_idx1,e_idx2)
    f_idx=np.intersect1d(f_idx1,f_idx2)

    e_and_f_idx = np.intersect1d(e_idx,f_idx)
    e_or_f_idx  = np.union1d(e_idx,f_idx)

    fig,ax = plt.subplots(4,1,figsize=(6,10),sharex=True)
    ax[0].plot(nsw+1,etot[nsw],'-',label='VASP')
    ax[0].plot(nsw+1,es[0],'-',label='nn')
    ax[1].plot(nsw+1,e_diff_per_atom,'-',label='VASP-NN')
    ax[2].plot(nsw+1, rmsd_f[nsw],'-',label='RMSE of force by nn and vasp')
    ax[0].set_ylabel('E (eV)');
    ax[1].set_ylabel(r'$\Delta E eV/atom$');
    ax[2].set_ylabel('RMSE force (eV/A)');

    ax[3].plot(e_or_f_idx+1,e_diff_per_atom[e_or_f_idx],'.',
             label='e ={0}-{1} || f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_or_f_idx)))
    ax[3].plot(e_and_f_idx+1,e_diff_per_atom[e_and_f_idx],'.',
             label='e ={0}-{1} && f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    ax[3].set_ylabel(r'$\Delta E eV/atom$');
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()

    np.savetxt(os.path.join(test_folder,'{0}_vid_e_or_f'.format(prefixs[0])),
               e_or_f_idx+1,fmt='%d',
               header='e ={0}-{1} || f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_or_f_idx)))
    np.savetxt(os.path.join(test_folder,'{0}_id_e_or_f'.format(prefixs[0])),
               e_or_f_idx,fmt='%d',
               header='e ={0}-{1} || f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_or_f_idx)))
    np.savetxt(os.path.join(test_folder,'{0}_vid_e_and_f'.format(prefixs[0])),
               e_and_f_idx+1,fmt='%d',
               header='e ={0}-{1} && f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    np.savetxt(os.path.join(test_folder,'{0}_id_e_and_f'.format(prefixs[0])),
               e_and_f_idx,fmt='%d',
               header='e ={0}-{1} && f = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    plt.savefig('vasp_vs_nn.png')
    fig.show()

elif mode == 'nn_only':
    max_dpgen_fs, mean_std_f = dev_nn(es,fs,vs,natoms)

    f_idx1 = np.where(abs(mean_std_f)<args.force_upper_cutoff)[0]
    f_idx2 = np.where(abs(mean_std_f)>args.force_lower_cutoff)[0]    
    f_idx=np.intersect1d(f_idx1,f_idx2)

    fig,ax = plt.subplots(3,1,figsize=(6,10),sharex=True)
#    ax[0].plot(nsw+1,etot[nsw],'-',label='VASP')
    for i in range(len(es)):
        ax[0].plot(es[i],'.',markersize=.5,alpha=.6,label=prefixs[i])
    dif_e = (es - nn[0])/natoms
    for i in range(len(es)):
        ax[1].plot(dif_e[i],'.',markersize=.5,alpha=.6,label=prefixs[i])
#    ax[2].plot(rmsd_f,'-',label='RMSE of force by nn and vasp)')
    ax[2].plot(mean_std_f,  label=' meand std f = {0}-{1}, {2}/{3}'.format(args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(f_idx), len(mean_std_f)))
    ax[2].plot([0,len(mean_std_f)],[args.force_upper_cutoff, args.force_upper_cutoff],'k--')
    ax[2].plot([0,len(mean_std_f)],[args.force_lower_cutoff, args.force_lower_cutoff],'k--')

    ax[0].set_ylabel('E (eV)');
    ax[1].set_ylabel(r'$\Delta E eV/atom$');
    ax[2].set_ylabel('model dev, force ave');
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    fig.savefig('nn_vs_nn.png')

    np.savetxt(os.path.join(test_folder,'model_dev_id'),
               f_idx,fmt='%d',
               header=' meand std f = {0}-{1}'.format(args.force_lower_cutoff, args.force_upper_cutoff))
    