#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020

same as std_vs_t.py, but output for different phases 

@author: jiedeng
"""

# extract energy 

# extract force

import dpdata.vasp.outcar
import argparse
import numpy as np
import os
from shared_functions import load_paths

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
    e_std  = np.std(energies[start:])
#    e_bar  = np.mean(energies)
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
    return e_std, f_std, v_std
    

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

def extract_sigma_vol_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    get_vol   = False
    get_sigma = False
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            get_sigma = True
        if 'volume of cell' in line:
            vol = float(line.split('volume of cell')[1].split(':')[1].split()[0])
            get_vol = True
        if get_vol and get_sigma:
            outcar_file.close()
            return vol, sigma


parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--outcar","-o",type=str,default = 'OUTCAR',help="name of outcar file")
parser.add_argument("--excudedfraction","-e",type=float,default = 0.2,help="excluded fraction, first few may not reach equil")
parser.add_argument("--suffix","-s",type=str,default = '',help="add recal or other suffix in the path")

args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,args.suffix) for path in tmp]

else:
    print("No folders point are provided. Use default value folders")
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

for path in paths:
    
    print(path)
    outcar = os.path.join(path,args.outcar)
#    outcar = '/Users/jiedeng/Documents/tmp/jd848/project_folder/ml/si-16/cont2/r3/OUTCAR'
    e_std,f_std,v_std = stds(outcar,args.excudedfraction)
    phase = get_phase(path)
#    print(e_std,f_std,v_std )
    vol, sigma = extract_sigma_vol_outcar(outcar)
#    print(vol, sigma)
    out.append([vol, sigma, e_std,f_std,v_std,phases[phase]])

out = np.array(out)
np.savetxt('std_vs_t_recal.out', out)


import matplotlib.pyplot as plt
plt.figure()
plt.plot(out[:,1],out[:,2],'o',label = 'e')
plt.plot(out[:,1],out[:,3],'^',label = 'f')
plt.plot(out[:,1],out[:,4],'s',label = 'v')
plt.xlabel('sigma (eV)')
plt.ylabel('standard deviation')
plt.yscale('log')
plt.legend()
plt.show()

