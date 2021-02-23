#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 12:46:26 2021

@author: jiedeng
"""

import pymatgen.io.vasp.outputs

filename = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/recal/OUTCAR'
out=pymatgen.io.vasp.outputs.Outcar(filename)

import pymatgen.io.xyz
xyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/recal/asap/ASAP-desc.xyz'
file=pymatgen.io.xyz.XYZ(xyz)

file=pymatgen.io.xyz.XYZ()


from pymatgen import Lattice, Structure, Molecule

file = Molecule.from_file(xyz)


struct.get_neighbor_list(r=2, sites = struct[:32]) # no detailed info which site corresponds to which
struct.get_neighbor_list(r=2, sites = struct[:2]) # no detailed info which site corresponds to which


struct=Structure.from_file('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/CONTCAR')

import pymatgen.analysis.local_env

brunnernn = pymatgen.analysis.local_env.BrunnerNN_real()
brunnernn.get_nn_info(struct,1)

