#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:58:13 2024

@author: jiedeng
"""

import argparse
from ase.io import read, write
from collections import Counter
# import numpy as np

def main(dump_file, output_file, element_order):
    # Read the trajectory
    traj = read(dump_file, index=':', format='lammps-dump-text')

    # Count the number of atoms of each type, assume atom type/number do not change
    type_counts = Counter(traj[0].get_atomic_numbers())
    # Generate the elements list based on user input and type counts
    elements_list = []
    for elem, count in zip(element_order, type_counts.values()):
        elements_list.extend([elem] * count)


    # Shift all atoms so the box starts at the origin
    for atoms in traj:
        cell_disp = atoms.get_celldisp()
        if len(atoms) == len(elements_list):
            atoms.set_chemical_symbols(elements_list)
            atoms.positions -= cell_disp
        else:
            raise ValueError("The length of elements_list does not match the number of atoms in the trajectory.")

    # Write the modified trajectory to an extended XYZ file
    write(output_file, traj, format='extxyz')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert LAMMPS dump to extended XYZ with specified elements.')
    parser.add_argument('--dump_file',default='npt.dump', help='Path to the LAMMPS dump file')
    parser.add_argument('--output_file',default='merge.xyz', help='Path for the output XYZ file')
    parser.add_argument('elements', nargs='+', help='List of element symbols in the order of their types in LAMMPS dump')
    
    args = parser.parse_args()
    main(args.dump_file, args.output_file, args.elements)
