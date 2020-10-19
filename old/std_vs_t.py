#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:41:34 2020
extract engery, forces
calculate std deviation of energy and force
plot


@author: jiedeng
"""

# extract energy 

# extract force

import dpdata.vasp.outcar
import argparse
import numpy as np
import os
from shared_functions import load_paths

def stds(outcar):
    """
    calculate std of outcar
    """
    atom_names, \
    atom_numbs, \
    atom_types, \
    cells, \
    coords, \
    energies, \
    forces, \
    tmp_virial = dpdata.vasp.outcar.get_frames(outcar)
    
    e_std  = np.std(energies)
#    e_bar  = np.mean(energies)
#    e_dif  = energies - e_bar
#    e_std2 = np.sqrt(sum(e_dif**2)/len(energies))
    
    f_bar = np.sum(forces,axis=0)/len(forces)   
    f_dif = forces - f_bar
    d_F   = np.sqrt((f_dif**2).sum(axis = 1).sum(axis = 1)/atom_numbs[0]/3) #eqn 13
    f_std = np.sqrt(sum(d_F**2)/len(d_F))

    v_bar = np.sum(tmp_virial,axis=0)/len(tmp_virial)   
    v_dif = tmp_virial - v_bar
    d_v   = np.sqrt((v_dif**2).sum(axis = 1).sum(axis = 1)/9) #eqn 13
    v_std = np.sqrt(sum(d_v**2)/len(d_v))
    
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


parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,'recal') for path in tmp]
#    paths = tmp
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd)
    paths = [cwd]

out = []
for path in paths:
    print(path)
    outcar = os.path.join(path,'outcar')
    e_std,f_std,v_std = stds(outcar)
    print(e_std,f_std,v_std )
    sigma = extract_sigma_outcar(outcar)
    print(sigma)
    out.append([sigma, e_std,f_std,v_std])

out = np.array(out)
np.savetxt('std_vs_t_recal.out', out)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(out[:,0],out[:,1],'o',label = 'e')
plt.plot(out[:,0],out[:,2],'^',label = 'f')
#plt.plot(out[:,0],out[:,3],'s',label = 'v')
plt.xlabel('sigma (eV)')
plt.ylabel('standard deviation')
plt.yscale('log')
plt.legend()
plt.show()

#import matplotlib.pyplot as plt
#plt.figure()
#plt.plot(out[:,0],out[:,1],'o',label = 'e')
#plt.plot(out[:,0],out[:,2],'^',label = 'f')
#plt.plot(out[:,0],out[:,3],'s',label = 'v')
#plt.legend()
