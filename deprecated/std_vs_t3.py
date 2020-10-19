#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020

same as std_vs_t.py, but output for different phases 

@author: jiedeng
"""

import dpdata.vasp.outcar
import argparse
import numpy as np
import os
from shared_functions import load_paths

from dp_test import test_ener, train_ener

import pandas as pd

def stds(outcar,exclude_fraction=0.2):
    """
    calculate std of outcar
    exclude_fraction : exclude the first 20% trajectory, not apply to recal case because the random choice of frame
    """
    atom_names, \
    atom_numbs, \
    atom_types, \
    cells, \
    coords, \
    energies, \
    forces, \
    tmp_virial = dpdata.vasp.outcar.get_frames(outcar)
    
    nsw   = len(energies)
    start = int(nsw*exclude_fraction)
    e_std = np.std(energies[start:])
    e_bar = np.mean(energies[start:])
#    e_dif  = energies - e_bar
#    e_std2 = np.sqrt(sum(e_dif**2)/len(energies))
    forces = forces[start:,]
    f_bar = np.sum(forces,axis=0)/len(forces)   
    f_dif = forces - f_bar
    d_F   = np.sqrt((f_dif**2).sum(axis = 1).sum(axis = 1)/atom_numbs[0]/3) #eqn 13
    f_std = np.sqrt(sum(d_F**2)/len(d_F))
    
    tmp_virial = tmp_virial[start:]
    v_bar = np.sum(tmp_virial,axis=0)/len(tmp_virial)   
    v_dif = tmp_virial - v_bar
    d_v   = np.sqrt((v_dif**2).sum(axis = 1).sum(axis = 1)/9) #eqn 13
    v_std = np.sqrt(sum(d_v**2)/len(d_v))
#    print(cell)
#    vol  = cell[0]
    return e_bar, v_bar, e_std, f_std, v_std, nsw
    

def extract_sigma_incar(incar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    with open(incar) as incar_file:
        for line in incar_file:
            if 'SIGMA' in line:
               sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
               return sigma

def extract_sigma_vol_kp_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    get_vol   = False
    get_sigma = False
    get_kinetic_pressure = False
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            get_sigma = True
        if 'volume of cell' in line:
            vol = float(line.split('volume of cell')[1].split(':')[1].split()[0])
            get_vol = True
        if 'kinetic pressure' in line:
            kp = float(line.split()[-2])
            get_kinetic_pressure = True
        if get_vol and get_sigma and get_kinetic_pressure:
            outcar_file.close()
            return vol, sigma, kp


parser = argparse.ArgumentParser()
### MUST SET###
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("-m", "--model", default="frozen_model.pb", type=str,help="Frozen model file to import")


parser.add_argument("--outcar","-o",type=str,default = 'OUTCAR',help="name of outcar file")
parser.add_argument("--excudedfraction","-e",type=float,default = 0.2,help="excluded fraction, first few may not reach equil")
parser.add_argument("--suffix","-s",type=str,default = '',help="add recal or other suffix in the path")

parser.add_argument("-S", "--set-prefix", default="set", type=str, 
                        help="The set prefix")
parser.add_argument("-n", "--numb-test", default=100, type=int, 
                        help="The number of data for test")
parser.add_argument("-r", "--rand-seed", type=int, 
                        help="The random seed")
parser.add_argument("--shuffle-test", action = 'store_true', 
                        help="Shuffle test data")
parser.add_argument("-d", "--detail-file", type=str, 
                        help="The file containing details of energy force and virial accuracy")


args   = parser.parse_args()
cwd    = os.getcwd()


print("**"*100)
print("**"*10,' '*10, 'ip format, MUST contain recal file' )
print("**"*10,' '*10, '/gpfs/loomis/project/kklee/jd848/ppv/10k/r1-519/recal/deepmd' )
print("**"*10,' '*10, '/gpfs/loomis/project/kklee/jd848/ppv/10k/r1-519/recal' )
print("**"*100)

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,args.suffix) for path in tmp]
else:
    print("No folders point are provided. Use current working path")
    inputpath = os.path.join(cwd)
    paths = [cwd]

out = []
phases = {'liquid':0,'pv':1,'ppv':2,'akimotoite':3,'enstatite':4, 'major':5, 'unkown':9}

def get_phase(path):
    if 'ppv' in path:
        phase = 'ppv'
    elif 'pv' in path:
        phase = 'pv'
    elif 'enstatite' in path:
        phase = 'enstatite'
    elif 'major' in path:
        phase = 'major'
    elif 'akimotoite' in path:
        phase = 'akimotoite'   
    else:
        print('phase unidentifiable')
        phase = 'unkown'
    return phase

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = tmp#[os.path.join(path,'recal') for path in tmp]
#    paths = tmp
else:
    print("No folders point are provided. Use current folder")
    inputpath = os.path.join(cwd)
    paths = [cwd]
    
inputs = {'model':args.model,
          'rand_seed':0,
          'system': None,
          'set_prefix':'set',
          'shuffle_test': args.shuffle_test,
          'numb_test':args.numb_test,
          'detail_file':args.detail_file}

eV_A3_2_GPa  = 160.21766208  # http://greif.geo.berkeley.edu/~driver/conversions.html
boltz        = 8.6173332e-5  

for path in paths:
    print('--'*40)
    print(path)
    print('--'*40)
    
    if  os.path.exists(os.path.join(path,'deepmd')):
        deepmd_path = os.path.join(path,'deepmd')
    elif os.path.exists(os.path.join(path,'deepmd_relax2')):
        deepmd_path = os.path.join(path,'deepmd_relax2')
    elif os.path.exists(os.path.join(path,'deepmd_relax')):
        deepmd_path = os.path.join(path,'deepmd_relax')
    else:
        deepmd_path = None
        
    if deepmd_path is None:
        pass
    else:
        
        recal_outcar = os.path.join(path,args.outcar)
        org_outcar   = os.path.normpath(os.path.join(path, os.pardir)) #os.path.join(path,'')
        org_outcar   = os.path.join(org_outcar,args.outcar)
    #    print(path)
        _, _, e_std_tr,f_std_tr,v_std_tr,train_size   = stds(recal_outcar)
        _, _, e_std_all,f_std_all,v_std_all,data_size = stds(org_outcar,0)
        e_bar, v_bar, e_std_eq,f_std_eq,v_std_eq,_    = stds(org_outcar,args.excudedfraction)
        
        
#        deepmd_path = os.path.join(path,'deepmd')
        
        inputs['system'] = deepmd_path
        
        num_test,sigma,natoms, l2e_test, l2ea_test, l2f_test, l2v_test = test_ener(inputs)
        num_tr,  sigma,natoms, l2e_tr,   l2ea_tr,   l2f_tr,   l2v_tr   = train_ener(inputs)
        
        phase = get_phase(path)
    #    print(e_std,f_std,v_std )
        vol, sigma, kp = extract_sigma_vol_kp_outcar(org_outcar)
        pii = np.mean(v_bar[:3])/vol*eV_A3_2_GPa
        p   = pii + kp/10  # in GPa
    #    print(vol, sigma)
        t   = round(sigma/boltz)
        l2v_test_gpa = l2v_test/vol*eV_A3_2_GPa
        l2v_tr_gpa   = l2v_tr/vol*eV_A3_2_GPa
    
        out.append([vol, t, p, e_bar, natoms, sigma, 
                    e_std_all, f_std_all, v_std_all,
                    e_std_eq, f_std_eq,v_std_eq,
                    e_std_tr, f_std_tr,v_std_tr,
                    l2ea_test, l2f_test, l2v_test,l2v_tr_gpa,
                    l2ea_tr, l2f_tr, l2v_tr,l2v_tr_gpa,
                    data_size, train_size,
                    phase, path])

columns = ['vol','T','P(GPa)','E(eV)','N','sigma(eV)',
           'Est_all','Fst_all','Vst_all',
           'Est_eq', 'Fst_eq', 'Vst_eq',
           'Est_tr', 'Fst_tr', 'Vst_tr',
           'Eaerr_test','Ferr_test','Verr_test','Verr_test(GPa)',
           'Eaerr_tr','Ferr_tr','Verr_tr','Verr_tr(GPa)',
           'Ndata','Ntrain',
           'phase','path']

out_pd = pd.DataFrame(out, columns = columns)

if not os.path.exists('out'):
    os.mkdir('out')
target = 'out'
out_pd.to_csv(os.path.join(target,'out_pd'))

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[1].plot(out_pd['T'],out_pd['P(GPa)'],'o')
ax[0].plot(out_pd['T'],out_pd['vol'],'o')
ax[2].plot(out_pd['T'],out_pd['E(eV)'],'o')
ax[0].set_xlabel('T'); ax[0].set_ylabel('P (GPa)')
ax[1].set_xlabel('T'); ax[1].set_ylabel('Vol (A3)')
ax[2].set_xlabel('T'); ax[2].set_ylabel('E (eV)')
fig.savefig(os.path.join(target,'pev.png'),bbox_inches='tight')

fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].scatter(out_pd['sigma(eV)'],out_pd['Est_eq'],c=out_pd['vol'])
ax[1].scatter(out_pd['sigma(eV)'],out_pd['Fst_eq'],c=out_pd['vol'])
f3 = ax[2].scatter(out_pd['sigma(eV)'],out_pd['Vst_eq'],c=out_pd['vol'])
plt.colorbar(f3)
ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Est (eV)')
ax[1].set_xlabel('sigma(eV)'); ax[1].set_ylabel('Fst (A3)')
ax[2].set_xlabel('sigma(eV)'); ax[2].set_ylabel('Vst (eV)')
fig.savefig(os.path.join(target,'simga_vs_std.png'),bbox_inches='tight')

fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].scatter(out_pd['sigma(eV)'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
ax[1].scatter(out_pd['sigma(eV)'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
f3=ax[2].scatter(out_pd['sigma(eV)'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
plt.colorbar(f3)
ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Eaerr/std')
ax[1].set_xlabel('sigma(eV)'); ax[1].set_ylabel('Fst/std')
ax[2].set_xlabel('sigma(eV)'); ax[2].set_ylabel('Vst/std')
fig.savefig(os.path.join(target,'simga_vs_err_std.png'),bbox_inches='tight')

fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].scatter(out_pd['Ndata'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
ax[1].scatter(out_pd['Ndata'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
f3=ax[2].scatter(out_pd['Ndata'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
plt.colorbar(f3)
ax[0].set_xlabel('Ndata'); ax[0].set_ylabel('Eaerr/std')
ax[1].set_xlabel('Ndata'); ax[1].set_ylabel('Fst/std')
ax[2].set_xlabel('Ndata'); ax[2].set_ylabel('Vst/std')
fig.savefig(os.path.join(target,'ndata_vs_err_std.png'),bbox_inches='tight')

