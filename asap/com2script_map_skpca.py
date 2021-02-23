#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 23:00:26 2020
This script shows how to translate the command line 

asap map -f per_atom.xyz -dm '[*]' -p map_per_atom --only_use_species 14 -ua skpca -n -1 -d 10

to script

Two ways to do that
1- relay on home-made map_tools, not recommended
2- relay on ASAP only


@author: jiedeng
"""

import numpy as np
from map_tools import map,  skpca

import matplotlib.pyplot as plt
from matplotlib import cm
from asaplib.data import ASAPXYZ
#from asaplib.reducedim import PCA
from asaplib.plot import *
import pandas as pd
import os


fxyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r1-k111/recal/asap/per_atom.xyz'
design_matrix = '[*]' 
prefix = 'ASPA'
output = 'chemiscope'
extra_properties = False
use_atomic_descriptors = True
only_use_species = 14
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
obj = {}        
obj1=map(obj, fxyz, design_matrix, prefix, output, extra_properties,
        use_atomic_descriptors, only_use_species, peratom, keepraw,
        color, color_column, color_label, colormap, color_from_zero, normalized_by_size,
        annotate, adjusttext, style, aspect_ratio)

scale =True  #False benchmark with the command line 
dimension = 10
axes = [0, 1]
kernel = 'linear'
kernel_parameter = None
sparse_mode = 'fps'
n_sparse = 100
n_sparse = -1
skpca(obj1,scale, dimension, axes,
          kernel, kernel_parameter, sparse_mode, n_sparse)



###########
########### The following is equivalent to above and have 
###########


asapxyz = ASAPXYZ(fxyz)
reduce_dict = {}
reduce_dict["preprocessing"] =  {"type": 'SCALE', 'parameter': None}
reduce_dict['skpca'] = {"type": 'SPARSE_KPCA',
                        'parameter':{"n_components": 10,
                                     "n_sparse": -1, # no sparsification
#                                     "scale":True,
                                "kernel": {"first_kernel": {"type": 'linear'}}}}

#reduce_dict['skpca'] = {"type": 'PCA',
#                        'parameter':{"n_components": 10,
#                         }}
from asaplib.reducedim import Dimension_Reducers
dreducer = Dimension_Reducers(reduce_dict)


dm  = asapxyz.get_atomic_descriptors(['SOAP-n6-l6-c6.0-g0.44'],14)
proj = dreducer.fit_transform(dm)

from asaplib.plot import Plotters

fig_spec = { 'outfile': None,
                'show': False,
                'title': None,
                'size': [8*1.1, 8],
                'cmap': 'gnuplot',
                    'components':{
                    'first_p': {'type': 'scatter', 'clabel': 'Relative enthalpy per TiO$_2$ [Kcal/mol]',
                                'vmin':None, 'vmax': None}
                    #'second_p': {"type": 'annotate', 'adtext': False} 
                    }
                   }
    
asap_plot = Plotters(fig_spec)
plotcolor = range(len(proj[:, [1,0]]))
asap_plot.plot(proj[:, [0,1]],plotcolor)

### plot command equivlaence
plt.scatter(proj[:, [0]],proj[:, [1]],c = plotcolor)


