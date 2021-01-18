#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 22:17:32 2020

@author: jiedeng
"""
import glob
from asaplib.data import ASAPXYZ
from ase.io import read, write
import numpy as np

class ASAPXYZs(ASAPXYZ):
    def __init__(self, fxyz=None, stride=1, periodic=True, fileformat=None):
        super(ASAPXYZ, self).__init__()
        if '*' in fxyz:
            self.fxyz = glob.glob(fxyz)
            print("Find matching input files with coordinates: ", self.fxyz)
        elif type(fxyz) is list:
            self.fxyz = fxyz
        else:
            self.fxyz = fxyz
        # essential
        self.stride = stride
        self.periodic = periodic
        if fileformat is not None:
            import ast
            self.fileformat = ast.literal_eval(fileformat)
        else:
            self.fileformat = {}

        # store the xyz file
        self.frames = None
        self.nframes = 0
        self.natom_list = []  # number of atoms for each frame
        self.total_natoms = 0  # total number of atoms for all frames
        self.global_species = []  # list of elements contains in all frames

        # record the state of the computation, e.g. which descriptors have been computed
        self.computed_desc_dict = {'data': {'fxyz': fxyz}}
        self.computed_desc_dict = {'descriptors': {}}
        # the conversion between tag of the descriptors and their acronyms
        self.tag_to_acronym = {'global': {}, 'atomic': {}}

        # we make a dictionary to store the computed descriptors
        self.global_desc = {}
        # this is for the atomic ones
        self.atomic_desc = {}

        # try to read the xyz file
        try:
            if isinstance(self.fxyz, (tuple, list)):
                self.frames = []
                for f in self.fxyz:
                    self.frames += read(f, slice(0, None, self.stride), **self.fileformat)
            else:
                self.frames = read(self.fxyz, slice(0, None, self.stride), **self.fileformat)
        except:
            raise ValueError('Exception occurred when loading the input file')

        self.nframes = len(self.frames)
        all_species = []
        for i, frame in enumerate(self.frames):
            # record the total number of atoms
            self.natom_list.append(len(frame.get_positions()))
            all_species.extend(frame.get_atomic_numbers())
            if not self.periodic or not np.sum(frame.get_cell()) > 0:
                frame.set_pbc([False, False, False])
            # we also initialize the descriptor dictionary for each frame
            self.global_desc[i] = {}
            self.atomic_desc[i] = {}

        self.total_natoms = np.sum(self.natom_list)
        self.max_atoms = max(self.natom_list)
        # Keep things in plain python for serialisation
        self.global_species = np.unique(all_species).tolist()
        print('load xyz file: ', self.fxyz,
              ', a total of ', str(self.nframes), 'frames',
              ', a total of ', str(self.total_natoms), 'atoms',
              ', with elements: ', self.global_species, '.')