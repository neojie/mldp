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
from map_tools import map,  skpca
    
obj = {}

import pandas as pd
#import shutil
import os
#out_pd = pd.read_excel('/Users/jiedeng/GD/ppt/2020/extreme_filter4_with_state.xlsx')
ppv = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'ppv')
pv = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'pv')
om8 = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'om8')
ak = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'ak')
en = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'en')
maj = pd.read_excel('/Users/jiedeng/GD/papers/pv_sum/dat_sum.xlsx',sheet_name = 'maj')

T0 =3000
om8_sel = om8[(om8['T'] == T0) & (om8['type'] =='equil') ]['local'].values
ppv_sel = ppv[(ppv['T'] == T0) & (ppv['P(GPa)']<300)]['local'].values
pv_sel  = pv[(pv['T'] == T0) & (pv['P(GPa)']<300)]['local'].values
ak_sel  = ak[(ak['T'] == T0) & (ak['P(GPa)']<300)]['local'].values
en_sel  = en[(en['T'] == T0) & (en['P(GPa)']<300)]['local'].values
maj_sel = maj[(maj['T'] == T0) & (maj['P(GPa)']<300)]['local'].values

#deepmds = out_pd['local_path'].values
file = 'ASAP-desc.xyz'
def get_vasp_asap_path(deepmds,file=file):
    fxyz_recal = []
    fxyz_vasp  = []
    for deepmd in deepmds:
        print(deepmd)
    #    if '3k' in deepmd:
        recal_dir  = os.path.abspath(os.path.join(deepmd, os.pardir))
        recal_asap = os.path.join(recal_dir,'asap')
        vaspdir    = os.path.abspath(os.path.join(recal_dir, os.pardir)) 
        vasp_asap  = os.path.join(vaspdir,'asap')
        if os.path.exists(os.path.join(recal_asap,file)):
            fxyz_recal.append(os.path.join(recal_asap,file))
        if os.path.exists(os.path.join(vasp_asap,file)):
            fxyz_vasp.append(os.path.join(vasp_asap,file))
    return fxyz_vasp, fxyz_recal
            



fxyz_vasp_om8, fxyz_recal_om8 = get_vasp_asap_path(om8_sel)
fxyz_vasp_ppv, fxyz_recal_ppv = get_vasp_asap_path(ppv_sel)
fxyz_vasp_pv, fxyz_recal_pv   = get_vasp_asap_path(pv_sel)
fxyz_vasp_ak, fxyz_recal_ak   = get_vasp_asap_path(ak_sel)
fxyz_vasp_en, fxyz_recal_en   = get_vasp_asap_path(en_sel)
fxyz_vasp_maj, fxyz_recal_maj = get_vasp_asap_path(maj_sel)

#fxyz = fxyz_vasp[:2]

fxyz = fxyz_recal_om8 +fxyz_recal_maj
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