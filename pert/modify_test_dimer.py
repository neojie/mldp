#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 22:13:57 2021

1- apply to dimer only
2-

/u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/extreme_filtered/re7/pv_gpu.pb
@author: jiedeng
"""

#from ase import Atoms
import ase.io.vasp
#import glob
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
import deepmd.DeepPot as DeepPot
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",default='0.vasp',help="input file")
parser.add_argument("--temperature","-t",type=float,default=4000,help="temperature")
parser.add_argument("--min_distance","-md",type=float,default=0.5,help="min distance")
parser.add_argument("--model","-m",type=str,default='/u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/mpmt/mm12/pv_Jan8_cpu.pb',help="model")

args   = parser.parse_args()

class DP_fparam(Calculator):
    name = "DP"
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, model, label="DP",fparam=0, **kwargs):
        Calculator.__init__(self, label=label, **kwargs)
        self.dp = DeepPot(model)
        self.type_dict = dict(zip(self.dp.get_type_map(), range(self.dp.get_ntypes())))
        self.fparam= [fparam]

    def calculate(self, atoms=None, properties=["energy", "forces", "stress"], system_changes=all_changes):
        coord = atoms.get_positions().reshape([1, -1])
        cell = atoms.get_cell().reshape([1, -1])
        symbols = atoms.get_chemical_symbols()
        atype = [self.type_dict[k] for k in symbols]
        e, f, v = self.dp.eval(coord, cell, atype,fparam=self.fparam)
        self.results['energy'] = e[0]
        self.results['forces'] = f[0]
        self.results['stress'] = v[0]

def cal_force_energy_dist(struct):
    force_vec = struct.get_forces()[0]
    force_mag = np.sqrt(np.sum(force_vec**2))
    energy = struct.get_potential_energy()
    dist = struct.get_all_distances()[0][1]
    return force_mag, energy, dist

kB = 8.617333e-5 # ev/k
fparam =  kB*args.temperature    
model      = args.model
calculator = DP_fparam(model=model,fparam=fparam)
struct    = ase.io.vasp.read_vasp(args.file)
struct.calc = calculator

pos0  = struct.positions[0] # atom 0 position
dist0 = struct.get_all_distances()[0][1]
abc = struct.get_cell_lengths_and_angles()[:3]
dists = np.linspace(args.min_distance,abc[0]/2,200)


force_mag = []
energy = []
distance = []
for dist in dists:
    pos1 = pos0+[dist,0,0]
    pos = [pos0,pos1]
    struct.set_positions(pos,True)
    _force_mag, _energy, _dist = cal_force_energy_dist(struct)
    force_mag.append(_force_mag)
    energy.append(_energy)
    distance.append(_dist)

import matplotlib.pyplot as plt
fig,ax = plt.subplots(2,1,figsize=(4,6),sharex=True)
ax[0].plot(distance,force_mag,'.')
ax[1].plot(distance,energy,'o')
ax[0].set_ylabel('force eV/A')
ax[1].set_ylabel('energy eV')
ax[1].set_xlabel('distance')
ax[0].grid()
ax[1].grid()
plt.show()