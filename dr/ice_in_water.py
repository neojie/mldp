#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 20:22:07 2020

@author: jiedeng
"""
#from asaplib.reducedim import KernelPCA





""" for maps and fits """

"""
/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/ASAP/build/lib/asaplib/cli/func_asap.py
"""
from bs_map import map,  skpca
    
obj = {}

#import pandas as pd
##import shutil
#import os
#out_pd = pd.read_excel('/Users/jiedeng/GD/ppt/2020/extreme_filter4_with_state.xlsx')
#deepmds = out_pd['local_path'].values
#
#file = 'ASAP-desc.xyz'
#fxyz_recal = []
#fxyz_vasp  = []
#for deepmd in deepmds:
#    if '3k' in deepmd:
#        recal_dir  = os.path.abspath(os.path.join(deepmd, os.pardir))
#        recal_asap = os.path.join(recal_dir,'asap')
#        vaspdir    = os.path.abspath(os.path.join(recal_dir, os.pardir)) 
#        vasp_asap  = os.path.join(vaspdir,'asap')
#        if os.path.exists(os.path.join(recal_asap,file)):
#            fxyz_recal.append(os.path.join(recal_asap,file))
#        if os.path.exists(os.path.join(vasp_asap,file)):
#            fxyz_vasp.append(os.path.join(vasp_asap,file))
#        

#fxyz = ['/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/recal/asap/ASAP-desc.xyz', 
#        '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/ASAP-desc.xyz']
#fxyz = ['/Users/jiedeng/GD/Learn/dscribe_learn/ice-in-water/ice-54/ice-54.xyz','/Users/jiedeng/GD/Learn/dscribe_learn/ice-in-water/liquid-1000/dataset_1000_eVAng.xyz']
fxyz = '/Users/jiedeng/GD/Learn/dscribe_learn/ice-in-water/notebooks/ASAP-desc.xyz'
design_matrix = '[*]'
prefix = 'ASPA'
output = 'chemiscope'
extra_properties = False
use_atomic_descriptors = False
only_use_species = None
peratom = False
keepraw = False
color = 'none'
color_column = 0
color_label = ''
colormap = 'gnuplot'
color_from_zero = False
normalized_by_size = True
annotate = 'none'
adjusttext = False
style = 'journal'
aspect_ratio =2
        
obj1=map(obj, fxyz, design_matrix, prefix, output, extra_properties,
        use_atomic_descriptors, only_use_species, peratom, keepraw,
        color, color_column, color_label, colormap, color_from_zero, normalized_by_size,
        annotate, adjusttext, style, aspect_ratio)


scale =True#False benchmark with the command line 
dimension = 5
axes = [0, 1]
kernel = 'linear'
kernel_parameter = None
sparse_mode = 'fps'
n_sparse = 100

skpca(obj1,scale, dimension, axes,
          kernel, kernel_parameter, sparse_mode, n_sparse)