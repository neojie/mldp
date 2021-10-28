#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 08:19:41 2021


# plot pressure, temperature, and energy without entropy (default)
# output the number of atoms

@author: jiedeng
"""
import os
import ase.io
import ase
import ase.io.vasp
import re
import numpy as np
import warnings
import argparse

# get nsw_tot 
parser = argparse.ArgumentParser(description="extract P E T V from OUTCAR")

parser.add_argument("--beg","-b", type=int,default=0, help="begin from index")
parser.add_argument("--end","-e", type=int,default=-1, help="end at index")
args = parser.parse_args()

def get_nsw_tot(vasp_dir = '.'):
    if os.path.exists(os.path.join(vasp_dir,'XDATCAR')):
        xdatcar = os.path.join(vasp_dir,'XDATCAR')
        for line in reversed(list(open(xdatcar))):
            tmp = line.rstrip()
            if "Direct" in tmp:
                nsw_tot = int(re.sub("[^0-9]", "",tmp))
                break
    elif os.path.exists(os.path.join(vasp_dir,'OSZICAR')):
        warnings.warn(
            "\n ****** XDATCAR does not exist ******\n"
            + "If nsw set > 9999, we cannot extract the results \n"
            + " And you should set nsw_tot manually!"
        )
        osizcar = os.path.join(vasp_dir,'OSZICAR')
        for line in reversed(list(open(osizcar))):
            tmp = line.rstrip()
            if tmp[0] ==' ' and "T=" in tmp:
                nsw_tot = int(tmp.split()[0])
                break
    # if OSIZCAR not exist, try OUTCAR
    elif os.path.exists(os.path.join(vasp_dir,'OUTCAR')):
        warnings.warn(
            "\n ****** XDATCAR does not exist ******\n"
            + "If nsw set > 9999, we cannot extract the results \n"
            + " And you should set nsw_tot manually!"
        )
        outcar = os.path.join(vasp_dir,'OUTCAR')
        for line in reversed(list(open(outcar))):
            tmp = line.rstrip()
            if "Iteration" in tmp:
                nsw_tot = int(tmp.split()[2].split('(')[0])
                break
    return nsw_tot

def read_outcar(vasp_dir,index):  
    outcar = os.path.join(vasp_dir,'OUTCAR')
    return ase.io.vasp.read_vasp_out(outcar,index)


nsw_tot = get_nsw_tot()
target_dir = '.'
vasp_dir   = '.'
outcar = os.path.join(vasp_dir,'OUTCAR')
atoms = []
dummy = ase.io.vasp.read_vasp_out(outcar)
natoms = dummy.get_global_number_of_atoms()

print("System:", dummy.get_chemical_formula(), dummy.get_global_number_of_atoms())


#property_list = []
#get_pos, get_energy = False, False
#nsw_i = 0
#ewithout = []

if 'pet.dat' in os.listdir():
    pass
else:
    from subprocess import call
    call("grep 'total pressure' OUTCAR | awk '{print $4}' > p.dat",shell=True)
    call("grep 'energy  without entropy' OUTCAR | awk '{print $4}' >e.dat",shell=True)
    call("grep '(temperature' OUTCAR | awk '{print $6}' > t.dat",shell=True)
    call("paste e.dat t.dat > temp",shell=True)
    call("paste p.dat temp >pet.dat",shell=True)
    call("rm e.dat t.dat p.dat",shell=True)

pet= np.loadtxt('pet.dat')
pet_be = pet[args.beg:args.end,:]
ave=np.mean(pet_be,axis=0)
print(f'{ave[0]/10}\t\t{ave[1]}\t{ave[2]}\t{dummy.get_volume()}')

import matplotlib.pyplot as plt
fig,ax = plt.subplots(3,1,figsize=(5,4),sharex=True)
ax[0].plot(pet[:,0]/10,label='P (GPa)')
ax[0].plot(args.beg, pet[args.beg,0]/10, 'ko')
ax[0].plot(nsw_tot, pet[args.end,0]/10, 'ko')
ax[1].plot(pet[:,1],label='E (eV)')
ax[2].plot(pet[:,2],label='T (K)')
ax[2].set_xlabel("Step")
ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()


#fig.savefig("stress.png")



#for i in range(len(lines)): # can start from a large number since typically OUTCAR first 5000 lines does not contatin that info
#    dummy=ase.io.vasp.read_vasp_out(outcar,i)
#    ewithout.append(dummy.get_potential_energy())  # energy without entropy


#for i in range(len(lines)): # can start from a large number since typically OUTCAR first 5000 lines does not contatin that info
#    line = lines[i]
#    if 'TOTAL-FORCE' in line:
#        pos_force = lines[i+2:i+atom_num+2]
#        pos_force = np.array([row.strip().split() for row in pos_force]).astype(np.float32)
#        pos    = pos_force[:,:3]
#        forces = pos_force[:,3:]
#        get_pos = True
#    if 'energy  without entropy=' in line:
#        # scalar must be stored as array instead
#        energy = np.array([float(line.split()[3])], dtype=np.float32) 
#        get_energy = True
#    
#    if get_pos and get_energy:
#        dummies[nsw_i].set_positions(pos)
#        property_list.append({'energy': energy,'forces':forces})
#        get_pos, get_energy = False, False
#        nsw_i += 1


"""
"""


#
#
#vasp_dir='.'
#get_nsw_tot(vasp_dir)
#out = read_outcar(vasp_dir,index)
#
#read_outcar(vasp_dir,1)
#
#ase.io.vasp.read_vasp_out(os.path.join(vasp_dir,'OUTCAR'),0)

#import pymatgen.io.vasp.outputs
#path = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv_hf_copy2/3k/solid1.5/r3-3k/'
#path = '/Users/jiedeng/Documents/tmp/jd848/Fe128Si12O24/melt4/r4-4000-FeSiO/'
##pymatgen.io.vasp.outputs.path(path+'/OUTCAR')
#pymatgen.io.vasp.outputs.Vasprun(path+'/vasprun.xml')
#1/(((2815-2623)+40)/(12*3500-250))



