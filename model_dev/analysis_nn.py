#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 10:48:49 2020
consider test and recal are in different folder
test resutls does not contain the orginal info
This is not a optimal way of practice. 
See analysis_recal.py for a better way

@author: jiedeng
"""

import numpy as np
import os
from model_dev_funcs import extract_lammps_log, extract_outcar, extract_nn_pred, dev_nn, dev_vasp_nn
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--test_folder","-tf",default='test',help="model test folder")
#parser.add_argument("--recal_foler","-rf",default='recal',help="recal folder")
parser.add_argument("--model_prefix","-mp",default='re4-re5-re6-re7-re8',help="recal folder")
#parser.add_argument("--energy_cutoff","-ec",default=0.0044,type=float,help="cutoff for energy")
parser.add_argument("--force_cutoff","-fc",default=0.27,type=float,help="cutoff for force")
parser.add_argument("--force_lower_cutoff","-flc",default=0.04,type=float,help="lower cutoff for force")
parser.add_argument("--force_upper_cutoff","-fuc",default=0.2,type=float,help="upper cutoff for force")
args   = parser.parse_args()
test_folder = args.test_folder
#recal_foler = args.recal_foler
prefixs = args.model_prefix.split('-')


### get confgis that were recalculated
#outcar = os.path.join(recal_foler,'OUTCAR')
natoms = 160
#etot, stress, forces, nsw = extract_outcar(outcar)
es,fs,vs                  = extract_nn_pred(test_folder,prefixs,natoms=160)
nn   = [np.mean(es,axis=0),np.mean(fs,axis=0),np.mean(vs,axis=0)]

#vasp = [etot,forces,stress]
#e_diff_per_atom = (etot[nsw] - nn[0])/160


#rmsd_e,rmsd_f,rmsd_v=dev_vasp_nn(vasp,nn,natoms=160)


if len(prefixs) > 1:
    max_dpgen_fs, mean_std_f = dev_nn(es,fs,vs)

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
#    ax[2].set_ylabel('RMSE force (eV/A)');
    ax[2].set_ylabel('model dev, force ave');
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    fig.savefig('corr.png')

    np.savetxt(os.path.join(test_folder,'model_dev_id'),
               f_idx,
               header=' meand std f = {0}-{1}'.format(args.force_lower_cutoff, args.force_upper_cutoff))


#    e_idx=np.where(abs(e_diff_per_atom)>0.0044)[0]
#    f_idx=np.where(abs(rmsd_f)>0.27)[0]
#    
#    cutoffs = np.linspace(0,max(abs(mean_std_f)))
#    ratios_e = []; num = []
#    ratios_f = [];  
#    for cutoff in cutoffs:
#        mod_dev_idx=np.where(abs(mean_std_f)>=cutoff)[0]
#        ratios_e.append(len(np.intersect1d(e_idx,mod_dev_idx))/len(mod_dev_idx))
#        ratios_f.append(len(np.intersect1d(f_idx,mod_dev_idx))/len(mod_dev_idx))
#        num.append(len(mod_dev_idx))
#    
#    fig1,ax1 = plt.subplots(1,1,figsize=(6,2),sharex=True)
#    ax1.plot(cutoffs,ratios_e,label='rmse e>nn error ')
#    ax1.plot(cutoffs,ratios_f,label='rmse f>nn error ')
#    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#    ax2.plot(cutoffs,num,'k')
#    ax1.legend()
#    ax1.set_xlabel('cutoff')
#    ax1.set_ylabel('correct frame fraction')
#    ax2.set_ylabel('tot frame selected')
#    ax1.grid()
#    fig1.savefig('cutoff.png')
#    
#    fig3,ax3 = plt.subplots(1,1,figsize=(6,2),sharex=True)
#    ax3.plot(f_idx,e_diff_per_atom[f_idx],'.',label='f > nn {0}'.format(len(f_idx)))
#    ef_idx = np.intersect1d(e_idx,f_idx)
#    ax3.plot(ef_idx,e_diff_per_atom[ef_idx],'.',label='f > nn && e > nn {0}'.format(len(ef_idx)))
#    ax3.legend()
#    fig3.savefig('idx_efc.png')
#    nsw_efc=np.array(nsw+1).astype(int)[ef_idx]
#    np.savetxt(os.path.join(recal_foler,'nsw_efc'),
#               nsw_efc,
#               header='energy cutoff {0} && force cutoff {1}'.format(args.energy_cutoff, args.force_cutoff))
#    fig3.show() 
#    fig.show()
#    fig1.show()
    

#else:
##    plt.figure()
##    plt.plot(rmsd_e/natoms)
##    plt.savefig('dif_e.png')
    


# find the optimum cutoff if we use threshold of 