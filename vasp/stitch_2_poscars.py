#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 12:43:02 2022

@author: jiedeng
"""
import numpy as np
import ase.io
from ase import Atoms 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--files","-f",type = str,nargs="+", default = ['POSCAR1',['POSCAR2']], help="files to be merged, default: POSCAR1 POSCAR2")
parser.add_argument("--gap1","-g1",type = float, default = 0.5, help="gap between the two poscar inside the slab, default: 0.5 A")
parser.add_argument("--gap2","-g2",type = float, default = 0.0, help=" the second gap on the edge of the slab considering PBC, default: 0.0 A")
parser.add_argument("--step","-s",type = int, default = 30, help="step, default: 30")
parser.add_argument("--axis","-a",type = int, default = 2, help="along which axis: default: 2 =>z; 0->x, 1->y")
parser.add_argument("--check_min_dis","-cmd", default = True, action='store_false', help="check min dis default:True")

args   = parser.parse_args()

file1,file2 = args.files[0],args.files[1]
gap =args.gap1
atom1 = ase.io.read(file1)
pos1 = atom1.positions
atom2 = ase.io.read(file2)
pos2 = atom2.positions

def getmins(atoms):
    dis = atoms.get_all_distances(mic=True)  ## mic consider the pbc
    
    mins = []
    for i in range(len(dis)-1):
        disi = dis[i]
        disi0 = disi[(i+1):]
        mins.append(min(disi0))
    return mins

def check_mins(atoms):
    mins=getmins(atoms)
    mindis = min(mins)
    if mindis<1:
        print("    !!! min distance = {0} < 1 A, consider increase the gap".format(mindis))
    else:
        print("    min distance is", mindis)

print("--- check min distance ----")
print(file1)
if args.check_min_dis:
    check_mins(atom1)
print("--- check min distance ----")
print(file2)
if args.check_min_dis:
    check_mins(atom2)
        


combined = atom1.get_chemical_symbols()+atom2.get_chemical_symbols() # note cannot just add the get_chemical_formula



axis = args.axis
pos2[:,axis] = pos2[:,axis]+gap+atom1.cell[axis][axis]
tmp=atom1.cell[:]
tmp[axis][axis] = tmp[axis][axis] + gap + atom2.cell[axis][axis] + args.gap2
    
newatom=Atoms(combined,positions=np.concatenate((pos1,pos2)),pbc=True,cell=tmp)


print("--- check min distance ----")
print('newposcar')
if args.check_min_dis:
    check_mins(newatom)
    
ase.io.write('newposcar.vasp',newatom,format='vasp')
# ase.io.write(atom1.get_chemical_formula(),atom1,format='vasp')
# new1atom1 = Atoms(atom1.get_chemical_symbols(),pos1,pbc=True,cell=atom1.cell)
# ase.io.write('newatom1',new1atom1,format='vasp')
