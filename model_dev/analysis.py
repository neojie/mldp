#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

analysis used for test without org dataset, analysis_recal for test with org dataset
Created on Tue Nov 24 10:48:49 2020


1) no outcar, >=1 dp test results, dp test results does not have original e, v, f
2) there is outcar, >=1 dp test results, dp test results does not have original e, v, f  => deprecated scenario
3) there is outcar, >=1 dp test results, dp test results does not have original e, v, f
See analysis_recal.py for a better way


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


1) no outcar, >=1 dp test results, dp test results does not have original e, v, f
2) there is outcar, >=1 dp test results, dp test results does not have original e, v, f  => deprecated scenario
3) there is outcar, >=1 dp test results, dp test results have original e, v, f

@author: jiedeng
"""

import numpy as np
import os
from model_dev_funcs import extract_org_nn_pred, dev_nn, dev_vasp_nn,extract_nn_pred
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="""
        example: analysis.py -tf recal -rf recal -mp mm4-mm5 
        extract average recal/mm4.e.out recal/mm5.e.out, the 1st col should be vasp values
        the 2nd col should be model predicted values. 
        mm4 and mm5 are averaged and compared against the vasp values.
        idx between given thresholds are extracted
        \n\n
        example: analysis.py -tf recal -mp mm4-mm5 
        extract recal/mm4.e.out recal/mm5.e.out, the 1st col is ignored
        the 2nd col should be model predicted values, mm4 and mm5 are caompared
        model devation is computed.
        idx between force thresholds are extracted
        \n\n
        Note these two modes, thresholds have different meanings
        """)

parser.add_argument("--test_folder","-tf",help="folder storing model test ")
parser.add_argument("--recal_foler","-rf",help="recal folder")
parser.add_argument("--model_prefix","-mp",default='mm2',help="recal folder, support one model or multiple models, e.g., re4-re5-re6")
parser.add_argument("--energy_lower_cutoff","-elc",default=0.0044*2,type=float,help="lower cutoff for energy")
parser.add_argument("--energy_upper_cutoff","-euc",default=0.1,type=float,help="upper cutoff for energy")
parser.add_argument("--force_lower_cutoff","-flc",default=0.27,type=float,help="lower cutoff for force")
parser.add_argument("--force_upper_cutoff","-fuc",default=1,type=float,help="upper cutoff for force")


args        = parser.parse_args()

#if not args.test_folder:
#    raise ValueError('No model tests are given')
mode = 'nn_only'
if args.recal_foler:
    recal_foler = args.recal_foler
    mode = 'nn_vasp'
    print("mode is {0}, nn average and vasp are compared".format(mode))
else:
    ## only nn test results are given
    test_folder = args.test_folder
    print("mode is {0}, nn model deviation is analyzed".format(mode))


#test_folder = args.test_folder
#recal_foler = args.recal_foler
prefixs     = args.model_prefix.split('-')


natoms = 160
if mode == 'nn_vasp':
    # nn and vasp results are stored togheter, no special nsw chosen needed
    es_org,fs_org,vs_org, es,fs,vs   = extract_org_nn_pred(test_folder,prefixs,natoms=natoms)
    etot,forces ,stress   = es_org[0],fs_org[0],vs[0]
    vasp = [etot,forces,stress]
    
    nsw = np.array(range(len(es_org[0])))
else:
    es,fs,vs               = extract_nn_pred(test_folder,prefixs,natoms=natoms)
    
nn   = [np.mean(es,axis=0),np.mean(fs,axis=0),np.mean(vs,axis=0)]


if mode == 'nn_vasp':
    e_diff_per_atom = (etot[nsw] - nn[0])/natoms
    rmsd_e,rmsd_f,rmsd_v=dev_vasp_nn(vasp,nn,natoms=natoms)

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
    plt.savefig('analysis_recal.png')
    fig.show()

elif mode == 'nn_only':
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
               f_idx,fmt='%d',
               header=' meand std f = {0}-{1}'.format(args.force_lower_cutoff, args.force_upper_cutoff))
    