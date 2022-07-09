#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 21:46:04 2022

@author: jiedeng
"""
import numpy as np
import ase.io
import MDAnalysis as mda
from stat_lib_base import cal_inter, fit_gds_double_sided,  temporal_mass_density_proj_coarse, pbc_x, unpbc_x, count_z

class Count(object):
    def __init__(self,
             # ase_xyz,
             # mda_xyz,
             xyz='merge.xyz',
             tagged_ele = 'He', # can be number or str
             mode='mass',
             nw=1, # for MgO-Fe system, we use 1, for Pv-H system, we use 2
             project_axis = 2,
             begin = 0,
             end = -1,
             step = 1,
             alpha = 2.5,
             verbose  = False,
             analyze = True):
        self.xyz = xyz
        self.ase_xyz = ase.io.read(xyz,index=':') 
        self.mda_xyz = mda.Universe(xyz)
        self.mode = mode
        self.nw = nw
        self.begin =begin
        if end == -1:
            self.end = len(self.ase_xyz)
        else:
            self.end = end
        self.step = step
        self.project_axis = project_axis
        self.alpha = alpha
        self.mesh = 1
        self.level = None
        self.verbose  = verbose
        
        self.num_atom = len(self.ase_xyz[0])
        
        self.ele_atomic_number = np.unique(self.ase_xyz[0].get_atomic_numbers())
        # ele_atomic_number = np.unique(tmp.get_atomic_numbers())

        self.ele_chemical_symbol = []
        for num in self.ele_atomic_number:
            self.ele_chemical_symbol.append(ase.data.chemical_symbols[num])
        print('System:', self.ase_xyz[0].get_chemical_formula())
        
        self.ele_indices = []
        for num in self.ele_atomic_number:
            self.ele_indices.append(np.array(range(self.num_atom))[self.ase_xyz[0].get_atomic_numbers()==num])
        # h_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==2]  #2
        # mg_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==12] # mg
        # si_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==26] # fe
        # o_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==8]  # o
        self.sum_elements_counts = []
        self.sum_phase_dimension = []
        self.sum_rechi =[]        
        if analyze:
            self.analyze()
            self.save_sum()