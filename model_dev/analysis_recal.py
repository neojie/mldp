#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

analysis used for test without org dataset, analysis_recal for test with org dataset
Created on Tue Nov 24 10:48:49 2020
1-Given NN0 where is trained with DS
2-MPMT  DataSet0=>DS0, DS0 has log file which contains potentail energy, and stress, dump file contains position (and atomic forces)
3-Use five or more different models to calculate model deviation, select idx0 to recal
  Select Rule : 1) mean of force model devi > 0.04, we do not impose upper band here
4-Recal using VASP for trajectories with idx0 to form VASP_DS0
  Note we have several pressure bins, so  VASP_DS0 = {VASP_DS01,VASP_DS02,VASP_DS03,...}
5-Compare energy, force, and virial of VASP_DS0 with those of NN0 prediction, select idx1 which is a subset of idx0
  Select Rule : 1)-force RMSE > reported error of NN && energy difference > reported error of NN
                2)- Skip several pressure bins, trained NN may be able to interpolate
6-Train with DS + VASP_DS0[idx1] to obtain NN1
7-Use NN1 to predict energy of whole VASP_DS0 dataset, select idx2
  Select Rule: 1)-force RMSE > reported error of NN && energy difference > reported error of NN
  in principle, idx2 should be a small submit of idx0, especially for those that were selected at step #5
8-Train with DS + VASP_DS0[idx1 || idx2] to obtain NN2
9-Use NN1 to predict energy of whole VASP_DS0 dataset, select idx3
  Select Rule: 1)-force RMSE > desired error of NN || energy differce > reported error of NN
10-Train with DS + VASP_DS0[idx1 || idx2 || idx3] to obtain NN3
reiterate step 8-9 until the # of outliars are very small amount, like 10% of total trajectory?

@author: jiedeng
"""

import numpy as np
import os
from model_dev_funcs import extract_org_nn_pred, dev_nn, dev_vasp_nn
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--test_folder","-tf",default='recal/deepmd',help="model test folder")
parser.add_argument("--recal_foler","-rf",default='recal',help="recal folder")
parser.add_argument("--model_prefix","-mp",default='re4-re5-re6-re7-re8',help="recal folder")
parser.add_argument("--energy_lower_cutoff","-elc",default=0.0044*2,type=float,help="lower cutoff for energy")
parser.add_argument("--energy_upper_cutoff","-euc",default=0.2,type=float,help="upper cutoff for energy")
parser.add_argument("--force_lower_cutoff","-fc",default=0.27,type=float,help="lower cutoff for force")
parser.add_argument("--force_upper_cutoff","-fuc",default=1,type=float,help="upper cutoff for force")

args        = parser.parse_args()
test_folder = args.test_folder
recal_foler = args.recal_foler
prefixs     = args.model_prefix.split('-')

#test_folder = '/Users/jiedeng/GD/papers/paperxx4_ml/sigma-20interval/40-60/test'
#recal_foler = '/Users/jiedeng/GD/papers/paperxx4_ml/sigma-20interval/40-60/recal'

### get confgis that were recalculated
outcar = os.path.join(recal_foler,'OUTCAR')
natoms = 160

es_org,fs_org,vs_org, es,fs,vs                  = extract_org_nn_pred(test_folder,prefixs,natoms=160)
etot,forces ,stress   = es_org[0],fs_org[0],vs[0]
nsw = range(len(es_org[0]))

nn   = [np.mean(es,axis=0),np.mean(fs,axis=0),np.mean(vs,axis=0)]
vasp = [etot,forces,stress]
e_diff_per_atom = (etot[nsw] - nn[0])/natoms


rmsd_e,rmsd_f,rmsd_v=dev_vasp_nn(vasp,nn,natoms=natoms)


if len(prefixs) > 1:
    # several models
    max_dpgen_fs, mean_std_f = dev_nn(es,fs,vs)

    dpgen_corr   = np.corrcoef(max_dpgen_fs,e_diff_per_atom)
    mean_corr    = np.corrcoef(mean_std_f,e_diff_per_atom)
    corr_this_f  = np.corrcoef(rmsd_f,mean_std_f)[0][1]
    corr_dpgen_f = np.corrcoef(rmsd_f,max_dpgen_fs)[0][1]

    fig,ax = plt.subplots(4,1,figsize=(6,10),sharex=True)
    ax[0].plot(nsw+1,etot[nsw],'-',label='VASP')
    ax[0].plot(nsw+1,es[0],'-',label='re4 nn')
    ax[1].plot(nsw+1,e_diff_per_atom,'-',label='VASP-NN')
    ax[2].plot(nsw+1, rmsd_f[nsw],'-',label='RMSE of force by nn and vasp)')
    ax[3].plot(nsw+1,mean_std_f[nsw],label=' corr(RMSE, model_dev): {:02.4f} / {:02.4f}'.format(corr_this_f,corr_dpgen_f))
    ax[0].set_ylabel('E (eV)');
    ax[1].set_ylabel(r'$\Delta E eV/atom$');
    ax[2].set_ylabel('RMSE force (eV/A)');
    ax[3].set_ylabel('model dev, force ave');
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    fig.savefig('corr.png')


    e_idx=np.where(abs(e_diff_per_atom)>0.0044)[0]
    f_idx=np.where(abs(rmsd_f)>0.27)[0]

    cutoffs = np.linspace(0,max(abs(mean_std_f)))
    ratios_e = []; num = []
    ratios_f = [];
    for cutoff in cutoffs:
        mod_dev_idx=np.where(abs(mean_std_f)>=cutoff)[0]
        ratios_e.append(len(np.intersect1d(e_idx,mod_dev_idx))/len(mod_dev_idx))
        ratios_f.append(len(np.intersect1d(f_idx,mod_dev_idx))/len(mod_dev_idx))
        num.append(len(mod_dev_idx))

    fig1,ax1 = plt.subplots(1,1,figsize=(6,2),sharex=True)
    ax1.plot(cutoffs,ratios_e,label='rmse e>nn error ')
    ax1.plot(cutoffs,ratios_f,label='rmse f>nn error ')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(cutoffs,num,'k')
    ax1.legend()
    ax1.set_xlabel('cutoff')
    ax1.set_ylabel('correct frame fraction')
    ax2.set_ylabel('tot frame selected')
    ax1.grid()
    fig1.savefig('cutoff.png')

    fig3,ax3 = plt.subplots(1,1,figsize=(6,2),sharex=True)
    ax3.plot(f_idx,e_diff_per_atom[f_idx],'.',label='f > nn {0}'.format(len(f_idx)))
    ef_idx = np.intersect1d(e_idx,f_idx)
    ax3.plot(ef_idx,e_diff_per_atom[ef_idx],'.',label='f > nn && e > nn {0}'.format(len(ef_idx)))
    ax3.legend()
    fig3.savefig('idx_efc.png')
    nsw_efc=np.array(nsw+1).astype(int)[ef_idx]
    np.savetxt(os.path.join(recal_foler,'nsw_efc'),
               nsw_efc,
               header='energy cutoff {0} && force cutoff {1}'.format(args.energy_cutoff, args.force_cutoff))
    fig3.show()
    fig.show()
    fig1.show()


else:
    # one model only

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
    ax[2].plot(nsw+1, rmsd_f[nsw],'-',label='RMSE of force by nn and vasp)')
    ax[0].set_ylabel('E (eV)');
    ax[1].set_ylabel(r'$\Delta E eV/atom$');
    ax[2].set_ylabel('RMSE force (eV/A)');

    ax[3].plot(e_or_f_idx+1,e_diff_per_atom[e_or_f_idx],'.',
             label='f ={0}-{1} || e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_or_f_idx)))
    ax[3].plot(e_and_f_idx+1,e_diff_per_atom[e_and_f_idx],'.',
             label='f ={0}-{1} && e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()

    nsw_efc=np.array(nsw+1).astype(int)[ef_idx]
    np.savetxt(os.path.join(test_folder,'vid_e_or_f'),
               e_or_f_idx+1,
               header='f ={0}-{1} && e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    np.savetxt(os.path.join(test_folder,'id_e_or_f'),
               e_or_f_idx,
               header='f ={0}-{1} && e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    np.savetxt(os.path.join(test_folder,'vid_e_and_f'),
               e_and_f_idx+1,
               header='f ={0}-{1} && e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    np.savetxt(os.path.join(test_folder,'id_e_and_f'),
               e_and_f_idx,
               header='f ={0}-{1} && e = {2}-{3}, {4}'.format(args.energy_lower_cutoff,args.energy_upper_cutoff,
                                                           args.force_lower_cutoff,args.force_upper_cutoff,
                                                           len(e_and_f_idx)))
    plt.savefig('analysis_recal.png')
    fig.show()
