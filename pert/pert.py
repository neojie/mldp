#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 13:01:27 2021

@author: jiedeng
"""

#import dpdata
import numpy as np
import ase.io.vasp
import copy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",help="input  file")
parser.add_argument("--num","-n",type=int,default=4,help="number of atoms to swap")

args   = parser.parse_args()


def swap(pos,idx0, idx1):
    pos_copy = copy.copy(pos)
    temp = copy.copy(pos_copy[idx0])
    pos_copy[idx0] = pos_copy[idx1]
    pos_copy[idx1] = temp
    return pos_copy


struct = ase.io.vasp.read_vasp(args.file)
pos = struct.positions
Zs = struct.get_atomic_numbers() 

n_mg = len(Zs[Zs==12]);  idx_mg = range(0,n_mg); 
n_si = len(Zs[Zs==14]);  idx_si = range(n_mg,n_mg+n_si)
n_o  = len(Zs[Zs==8]);   idx_o  = range(n_mg+n_si,n_mg+n_si+n_o)

num_swap = args.num; 
idx_mg_mgo_swap    =  np.random.choice(idx_mg, num_swap, replace=False)
idx_o_mgo_swap     =  np.random.choice(idx_o, num_swap, replace=False)

num_swap = args.num; 
idx_si_sio_swap    =  np.random.choice(idx_si, num_swap, replace=False)
idx_o_sio_swap     =  np.random.choice(idx_o, num_swap, replace=False)

num_swap = args.num; 
idx_mg_mgsi_swap    =  np.random.choice(idx_mg, num_swap, replace=False)
idx_si_mgsi_swap    =  np.random.choice(idx_si, num_swap, replace=False)


pos_mgo = swap(pos,idx_mg_mgo_swap,idx_o_mgo_swap)
pos_sio = swap(pos,idx_si_sio_swap,idx_o_sio_swap)
pos_mgsi = swap(pos,idx_mg_mgsi_swap,idx_si_mgsi_swap)

import dpdata 

def save(pos,vasp='pos_mgo.vasp',lmp='pos_mgo.lmp'):
    struct.set_positions(pos)  
    ase.io.write(vasp,struct,format = 'vasp',vasp5=True,direct=True)    
    ls=dpdata.System(vasp,fmt='vasp/poscar')
    ls.to_lammps_lmp(lmp)    

save(pos_mgo,'pos_mgo.vasp','pos_mgo.lmp')
save(pos_sio,'pos_sio.vasp','pos_sio.lmp')
save(pos_mgsi,'pos_mgsi.vasp','pos_mgsi.lmp')


