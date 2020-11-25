#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020
shared function, select what to drop
1-analysis deepmd only, analysis outcar only
outcar only, if no model provided

@author: jiedeng
"""

print("**"*40)
print("*"*15,' '*19, 'ip format',' '*15, "*"*18 )
print("*"*15,' '*10, 'xx, xx/recal/xx/recal/deepmd*',' '*10,"*"*12 )
print("*"*15,' '*10, '  if no ip, in recal folder  ',' '*10,"*"*12 )

print("**"*40)

import dpdata.vasp.outcar
import argparse
import numpy as np
import os
from shared_functions import load_paths

import pandas as pd
import glob

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
    tmp_virial = dpdata.vasp.outcar.get_frames(outcar)   ## Note here tmp_varial is in kB!!!
    
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
    f_std = np.sqrt(sum(d_F**2)/len(d_F)) # => should be len(d_F) -1
    
    tmp_virial = tmp_virial[start:]
    vol = (np.cross(cells[0][0], cells[0][1])*cells[0][2]).sum()
#    print("cells", cells)
#    print("vol is",vol)
    tmp_virial = tmp_virial/10*vol/eV_A3_2_GPa   # kB to eV
    
    v_bar = np.sum(tmp_virial,axis=0)/len(tmp_virial)   
    v_dif = tmp_virial - v_bar
    d_v   = np.sqrt((v_dif**2).sum(axis = 1).sum(axis = 1)/9) #eqn 13
    v_std = np.sqrt(sum(d_v**2)/len(d_v))
#    print(cell)
#    vol  = cell[0]
    return e_bar, v_bar, e_std, f_std, v_std, nsw,sum(atom_numbs)


def stds2(outcar,exclude_fraction=0.2):
    """
    calculate std of outcar
    exclude_fraction : exclude the first 20% trajectory, not apply to recal case because the random choice of frame
    for some case, we exclude some set, so the nsw should be smaller than that in the OUTCAR
    """
    atom_names, \
    atom_numbs, \
    atom_types, \
    cells, \
    coords, \
    energies, \
    forces, \
    tmp_virial = dpdata.vasp.outcar.get_frames(outcar)   ## Note here tmp_varial is in kB!!!
    
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
    vol = (np.cross(cells[0][0], cells[0][1])*cells[0][2]).sum()
#    print("cells", cells)
#    print("vol is",vol)
    tmp_virial = tmp_virial/10*vol/eV_A3_2_GPa   # kB to eV
    
    v_bar = np.sum(tmp_virial,axis=0)/len(tmp_virial)   
    v_dif = tmp_virial - v_bar
    d_v   = np.sqrt((v_dif**2).sum(axis = 1).sum(axis = 1)/9) #eqn 13
    v_std = np.sqrt(sum(d_v**2)/len(d_v))
#    print(cell)
#    vol  = cell[0]
    
    ### Ntrain
    energies=glob.glob(deepmd_path+'/set*'+'/energy.npy')
    count = 0
    for energy in energies:
        count += len(np.load(energy))
    return e_bar, v_bar, e_std, f_std, v_std, count,sum(atom_numbs)
    
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
parser.add_argument("-m", "--model", type=str,help="Frozen model file to import, no model specify => no model mode")
parser.add_argument("-mo", "--mo", type=str,help="model only")
parser.add_argument("-de", "--deepmd", type=str,default = 'deepmd',help="deepmd folder, could be deepmd-deepmd_relax, use - as separator")


parser.add_argument("--outcar","-o",type=str,default = 'OUTCAR',help="name of outcar file")
parser.add_argument("--excudedfraction","-e",type=float,default = 0.2,help="excluded fraction, first few may not reach equil")
#parser.add_argument("--suffix","-s",type=str,default = '',help="add recal or other suffix in the path")

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

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath,level='recal')
#    paths = [os.path.join(path,args.suffix) for path in tmp]
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


if args.model:  ##
    from dp_test import test_ener, train_ener
    inputs = {'model':args.model,
              'rand_seed':0,
              'system': None,
              'set_prefix':'set',
              'shuffle_test': args.shuffle_test,
              'numb_test':args.numb_test,
              'detail_file':args.detail_file}

eV_A3_2_GPa  = 160.21766208  # http://greif.geo.berkeley.edu/~driver/conversions.html
boltz        = 8.6173332e-5  

import random
tag = str(random.randint(0,1e5))  ## put a tag to avoid data overprint
for path in paths:
    print('--'*40)
    print(path)
    print('--'*40)
    ## check the source data, some only has deepmd, not deepmd_relax2
    ## and vice versa. deepmd_relax2 is for relax process only
    ## to do the statisitics, we only care where there is deepmd
    if  os.path.exists(os.path.join(path,args.deepmd)): 
        deepmd_path = os.path.join(path,args.deepmd)     
        recal_outcar = os.path.join(path,args.outcar)
        org_outcar_tmp   = os.path.normpath(os.path.join(path, os.pardir)) #os.path.join(path,'')
        org_outcar       = os.path.join(org_outcar_tmp,args.outcar)
    #    print(path)
        _, _, e_std_tr,f_std_tr,v_std_tr,train_size,_   = stds2(recal_outcar)
        _, _, e_std_all,f_std_all,v_std_all,data_size,_ = stds(org_outcar,0)
        e_bar, v_bar, e_std_eq,f_std_eq,v_std_eq,_,natoms    = stds(org_outcar,args.excudedfraction)   
        
        phase = get_phase(path)
    #    print(e_std,f_std,v_std )
        vol, sigma, kp = extract_sigma_vol_kp_outcar(org_outcar)
        pii = (v_bar[0][0] + v_bar[1][1] + v_bar[2][2])/3/vol*eV_A3_2_GPa
        p   = pii + kp/10  # in GPa
        t   = round(sigma/boltz)
#        print("v_bar is:",v_bar)
#        print("kp is:",kp)
        if args.model:                
            inputs['system'] = deepmd_path        
            num_test,sigma,natoms, l2e_test, l2ea_test, l2f_test, l2v_test = test_ener(inputs)
            num_tr,  sigma,natoms, l2e_tr,   l2ea_tr,   l2f_tr,   l2v_tr   = train_ener(inputs)
        
        #    print(vol, sigma)
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
        else:
            out.append([vol, t, p, e_bar, natoms, sigma, 
                        e_std_all, f_std_all, v_std_all,
                        e_std_eq, f_std_eq,v_std_eq,
                        e_std_tr, f_std_tr,v_std_tr,
                        data_size,
                        phase, path])            

    
    
if args.model:   
    columns = ['vol','T','P(GPa)','E(eV)','N','sigma(eV)',
               'Est_all','Fst_all','Vst_all',
               'Est_eq', 'Fst_eq', 'Vst_eq',
               'Est_tr', 'Fst_tr', 'Vst_tr',
               'Eaerr_test','Ferr_test','Verr_test','Verr_test(GPa)',
               'Eaerr_tr','Ferr_tr','Verr_tr','Verr_tr(GPa)',
               'Ndata','Ntrain',
               'phase','path']
else:
    columns = ['vol','T','P(GPa)','E(eV)','N','sigma(eV)',
               'Est_all','Fst_all','Vst_all',
               'Est_eq', 'Fst_eq', 'Vst_eq',
               'Est_tr', 'Fst_tr', 'Vst_tr',
               'Ndata',
               'phase','path']    

out_pd = pd.DataFrame(out, columns = columns)

if not os.path.exists('out'):
    os.mkdir('out')
target = 'out'
if args.model:
    out_pd.to_excel(os.path.join(target,'out_pd_model'+'_'+tag+'.xlsx'))
else:
    out_pd.to_excel(os.path.join(target,'no_model'+'_'+tag+'.xlsx'))
    
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].plot(out_pd['T'],out_pd['P(GPa)'],'o')
ax[1].plot(out_pd['T'],out_pd['vol'],'o')
ax[2].plot(out_pd['T'],out_pd['E(eV)'],'o')
ax[0].set_xlabel('T'); ax[0].set_ylabel('P (GPa)')
ax[1].set_xlabel('T'); ax[1].set_ylabel('Vol (A3)')
ax[2].set_xlabel('T'); ax[2].set_ylabel('E (eV)')
fig.savefig(os.path.join(target,'pev_{0}.png'.format(tag)),bbox_inches='tight')

fig, ax = plt.subplots(1,3,figsize=(14,4))
ax[0].scatter(out_pd['sigma(eV)'],out_pd['Est_eq'],c=out_pd['vol'])
ax[1].scatter(out_pd['sigma(eV)'],out_pd['Fst_eq'],c=out_pd['vol'])
f3 = ax[2].scatter(out_pd['sigma(eV)'],out_pd['Vst_eq'],c=out_pd['vol'])
plt.colorbar(f3)
ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Est (eV/atom)')
ax[1].set_xlabel('sigma(eV)'); ax[1].set_ylabel('Fst (eV/A)')
ax[2].set_xlabel('sigma(eV)'); ax[2].set_ylabel('Vst (eV/atom)')
fig.savefig(os.path.join(target,'simga_vs_std_{0}.png'.format(tag)),bbox_inches='tight')

if args.model:
    fig, ax = plt.subplots(1,3,figsize=(14,4))
    ax[0].scatter(out_pd['sigma(eV)'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
    ax[1].scatter(out_pd['sigma(eV)'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
    f3=ax[2].scatter(out_pd['sigma(eV)'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
    plt.colorbar(f3)
    ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Eaerr/std')
    ax[1].set_xlabel('sigma(eV)'); ax[1].set_ylabel('Fst/std')
    ax[2].set_xlabel('sigma(eV)'); ax[2].set_ylabel('Vst/std')
    fig.savefig(os.path.join(target,'simga_vs_err_std_{0}.png'.format(tag)),bbox_inches='tight')
    
    fig, ax = plt.subplots(1,3,figsize=(14,4))
    ax[0].scatter(out_pd['Ndata'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
    ax[1].scatter(out_pd['Ndata'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
    f3=ax[2].scatter(out_pd['Ndata'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
    plt.colorbar(f3)
    ax[0].set_xlabel('Ndata'); ax[0].set_ylabel('Eaerr/std')
    ax[1].set_xlabel('Ndata'); ax[1].set_ylabel('Fst/std')
    ax[2].set_xlabel('Ndata'); ax[2].set_ylabel('Vst/std')
    fig.savefig(os.path.join(target,'ndata_vs_err_std_{0}.png'.format(tag)),bbox_inches='tight')

