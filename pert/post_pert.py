#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 00:21:17 2021

@author: jiedeng
"""

import ase.io.lammpsdata
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",default='mgsio3.dump',help="input  file")
parser.add_argument("--format","-ft",default='dump',help="input  file")

args   = parser.parse_args()

file = args.file
#file = '/Users/jiedeng/GD/papers/paperxx4_ml/post_nn/generate_new_poscar/new/re_4k_p2/mgsio3.dump'

cutoffs = {'MgMg': 1.2, 'MgSi': 1.15, 'MgO': 1.1, 
           'SiSi': 1.2, 'SiO':  1.1, 'OO': 0.9}

if args.format == 'dump':
    struct = ase.io.read(file,format='lammps-dump-text',index=0)
elif args.format == 'vasp':
    struct = ase.io.read(file,format='vasp')
elif args.format == 'lmp':
    struct = ase.io.read(file,format='lammps-data',style='atomic')
else:
    print('format not supported!')
  
pos = struct.positions
Zs = struct.get_atomic_numbers() 

n_mg = len(Zs[Zs==Zs[0]]);  idx_mg = range(0,n_mg); 
n_si = len(Zs[Zs==Zs[1]]);  idx_si = range(n_mg,n_mg+n_si)
n_o  = len(Zs[Zs==Zs[2]]);   idx_o  = range(n_mg+n_si,n_mg+n_si+n_o)

def get_min_dist(dist,idx0,idx1):
    mindist = np.inf
    for i in idx0:
        for j in idx1:
            if i < j and dist[i,j] <mindist:
#                print(dist[i,j])
                mindist = dist[i,j]   
    print('***'*10)
    print(i,j,mindist)
    return mindist

#struct = ase.io.read('/Users/jiedeng/GD/papers/paperxx4_ml/post_nn/generate_new_poscar/50000.vasp',format='vasp')
#get_min_dist(struct.get_all_distances(mic=True),range(160),range(160))

idx = []
dist = struct.get_all_distances(mic=True)
for ii in range(100000):
    try:
        print('--'*20)
        print(ii)
        
        if args.format == 'dump':
            struct = ase.io.read(file,format='lammps-dump-text',index=ii)
        elif args.format == 'vasp':
            struct = ase.io.read(file,format='vasp')
        elif args.format == 'lmp':
            struct = ase.io.read(file,format='lammps-data',style='atomic')
        else:
            print('format not supported!')

        dist=struct.get_all_distances(mic=True)
        mindist = get_min_dist(dist,idx_mg,idx_mg)             
        if mindist < cutoffs['MgMg']:
            print('MgMg ', mindist)
            raise ValueError('MgMg ', mindist)

        mindist = get_min_dist(dist,idx_mg,idx_si)             
        if mindist < cutoffs['MgSi']:
            print('MgSi ', mindist)
            raise ValueError('MgSi ', mindist)

        mindist = get_min_dist(dist,idx_mg,idx_o)             
        if mindist < cutoffs['MgO']:
            print('MgO ', mindist)
            raise ValueError('MgO ', mindist)

        mindist = get_min_dist(dist,idx_si,idx_si)             
        if mindist < cutoffs['SiSi']:
            print("SiSi ", mindist)
            raise ValueError("SiSi ", mindist)

        mindist = get_min_dist(dist,idx_si,idx_o)             
        if mindist < cutoffs['SiO']:
            print("SiO ", mindist)
            raise ValueError("SiO ", mindist)

        mindist = get_min_dist(dist,idx_o,idx_o)             
        if mindist < cutoffs['OO']:
            print("OO ", mindist)
            raise ValueError("OO ", mindist)
            
        idx.append(ii)
    except:
        print(ii)
        break
        
print('done')
print(idx)