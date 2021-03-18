#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 20:47:45 2021

https://pymatgen.org/pymatgen.io.vasp.outputs.html

@author: jiedeng
"""

from ase import neighborlist
from ase.build import molecule
from scipy import sparse
mol = molecule('CH3CH2OH')
cutOff = neighborlist.natural_cutoffs(mol)
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
neighborList.update(mol)
matrix = neighborList.get_connectivity_matrix()
#or: matrix = neighborlist.get_connectivity_matrix(neighborList.nl)
n_components, component_list = sparse.csgraph.connected_components(matrix)
idx = 1
molIdx = component_list[idx]
print("There are {} molecules in the system".format(n_components))
print("Atom {} is part of molecule {}".format(idx, molIdx))
molIdxs = [ i for i in range(len(component_list)) if component_list[i] == molIdx ]
print("The following atoms are part of molecule {}: {}".format(molIdx, molIdxs))




import io
from ase.io import read
def xyz2ase(xyz_str):
    """    Convert a xyz file to an ASE atoms object via in-memory file (StringIO).    """
    xyzfile = io.StringIO()
    xyzfile.write(xyz_str)
    mol = read(xyzfile, format="xyz")
    return mol

xyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/ASAP-desc.xyz'
at=xyz2ase(xyz)

import ase.io.vasp 
out =  ase.io.vasp.read_vasp_out('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/OUTCAR',index=0)

from ase.neighborlist import NeighborList
nl = NeighborList(2*[5])
nl.update(out)
indices, offsets = nl.get_neighbors(0)