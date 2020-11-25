#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 10:48:49 2020

@author: jiedeng
"""

import numpy as np
import os
from model_dev_funcs import extract_lammps_log, extract_outcar, extract_nn_pred, dev_nn, dev_vasp_nn
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--logfile","-log",default='rerun.lammps.output',help="lammps log file")
parser.add_argument("--test_folder","-tf",default='test',help="model test folder")
parser.add_argument("--recal_foler","-rf",default='recal',help="recal folder")
parser.add_argument("--skip_header","-sh",default=125,type=int,help="# of lines skipped of logfile header")
parser.add_argument("--skip_footer","-sf",default=16,type=int,help="# of lines skipped of logfile footer")


args   = parser.parse_args()
logfile     = args.logfile
test_folder = args.test_folder
recal_foler = args.recal_foler
skip_header = args.skip_header
skip_footer = args.skip_footer


logfile     = '/Users/jiedeng/GD/papers/paperxx4_ml/sigma-20interval/40-60/rerun.lammps.output'
test_folder = '/Users/jiedeng/GD/papers/paperxx4_ml/sigma-20interval/40-60/test'
recal_foler = '/Users/jiedeng/GD/papers/paperxx4_ml/sigma-20interval/40-60/recal'

### get confgis that were recalculated
outcar = os.path.join(recal_foler,'OUTCAR')

etot, stress, forces, nsw = extract_outcar(outcar)
#out                       = extract_lammps_log(logfile,skip_header,skip_footer)# the pb was old version
es,fs,vs                  = extract_nn_pred(test_folder,['re4','re5','re6','re7','re8'],natoms=160)

e_diff_per_atom = (etot[nsw] - es[0])/160

max_dpgen_fs, mean_std_f = dev_nn(es,fs,vs)

dpgen_corr = np.corrcoef(max_dpgen_fs,e_diff_per_atom)
mean_corr  = np.corrcoef(mean_std_f,e_diff_per_atom)


fig,ax = plt.subplots(4,1,figsize=(6,10),sharex=True)
ax[0].plot(nsw+1,etot[nsw],'-',label='VASP')
#ax[0].plot(nsw+1,out[nsw,2],'-',label='NN')
ax[0].plot(nsw+1,es[0],'-',label='re4 nn')
ax[1].plot(nsw+1,e_diff_per_atom,'-',label='VASP-NN')
ax[2].plot(nsw+1,max_dpgen_fs,label='dpgen, corr = {:02.4f}'.format(dpgen_corr[0][1]))
ax[3].plot(nsw+1,mean_std_f,label='mean of std, corr =  {:02.4f}'.format(mean_corr[0][1]))
ax[0].set_ylabel('E (eV)');
ax[1].set_ylabel(r'$\Delta E eV/atom$');
ax[2].set_ylabel('model dev, force max');
ax[3].set_ylabel('model dev, force ave');
ax[0].legend()
ax[1].legend()
ax[2].legend()
ax[3].legend()
fig.savefig('corr.png')
fig.show()

fig,ax = plt.subplots(4,1,figsize=(6,10),sharex=True)
ax[0].plot(nsw+1,etot[nsw],'-',label='VASP')
#ax[0].plot(nsw+1,out[nsw,2],'-',label='NN')
ax[0].plot(nsw+1,es[0],'-',label='re4 nn')
ax[1].plot(nsw+1,e_diff_per_atom,'-',label='VASP-NN')
ax[2].plot(nsw+1,max_dpgen_fs,label='dpgen, corr = {:02.4f}'.format(dpgen_corr[0][1]))
ax[3].plot(nsw+1,mean_std_f,label='mean of std, corr =  {:02.4f}'.format(mean_corr[0][1]))
ax[0].set_ylabel('E (eV)');
ax[1].set_ylabel(r'$\Delta E eV/atom$');
ax[2].set_ylabel('model dev, force max');
ax[3].set_ylabel('model dev, force ave');
ax[0].legend()
ax[1].legend()
ax[2].legend()
ax[3].legend()
fig.savefig('corr.png')
fig.show()


vasp = [etot,forces,stress]
nn   = [es[0],fs[0],vs[0]]

rmsd_e,rmsd_f,rmsd_v=dev_vasp_nn(vasp,nn,natoms=160)


fig,ax = plt.subplots(2,1,figsize=(6,5),sharex=True)
corr_this_f  = np.corrcoef(rmsd_f,mean_std_f)[0][1]
corr_dpgen_f = np.corrcoef(rmsd_f,max_dpgen_fs)[0][1]
corr_this_v  = np.corrcoef(rmsd_v,mean_std_f)[0][1]
corr_dpgen_v = np.corrcoef(rmsd_v,max_dpgen_fs)[0][1]
ax[0].plot(nsw+1, rmsd_f[nsw],'-',label='RMSE corr: {:02.4f} / {:02.4f}'.format(corr_this_f,corr_dpgen_f))
ax[1].plot(nsw+1, rmsd_v[nsw],'-',label='RMSE corr: {:02.4f} / {:02.4f}'.format(corr_this_v,corr_dpgen_v))
ax[0].legend()
ax[1].legend()
ax[0].set_ylabel('RMSE force (eV/A)');
ax[1].set_ylabel('RMSE virial (eV/A)');



var_e = np.var(es,axis=0);std_e = np.std(es,axis=0)
var_f = np.var(fs,axis=0);std_f = np.std(fs,axis=0) 
max_var_f = np.max(var_f,axis=1);max_std_f = np.max(std_f,axis=1)
var_v = np.var(vs,axis=0);std_v = np.std(vs,axis=0);
max_var_v = np.max(var_v,axis=1);max_std_v = np.max(std_v,axis=1)


np.corrcoef(var_e,e_diff_per_atom)
np.corrcoef(max_var_v,e_diff_per_atom)

np.corrcoef(max_var_v,e_diff_per_atom)
plt.plot(max_var_v)
plt.plot(e_diff_per_atom)
#
#plt.plot((etot - es[0])/160)
#
#len(np.where(abs((etot - es[0])/160) >0.01)[0])
#len(np.where(abs((etot - out[nsw,2])/160) >0.004)[0])
#
#plt.plot(es[0],out[nsw,2])
#
#
#np.corrcoef(rmsd_f,mean_std_f)
