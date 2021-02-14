#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 23:00:26 2020

Refs 
com2script_map_skpca 
/Users/jiedeng/GD/Learn/dscribe_learn/Mapping-the-space-of-materials-and-molecules/TiO2/TiO2-ASAP-KPCA.py

@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
from asaplib.data import ASAPXYZ
import pandas as pd
import os

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
file = 'per_atom.xyz'
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


#get total frames from corresponding oucar file for recal folder
def get_nframes(deepmds,dir_type = 'deepmd'):
    nsws = []
    for deepmd in deepmds: # deepmds can be asap
        print(deepmd)
        if dir_type == 'deepmd' or dir_type == 'asap':
            recal_dir = os.path.abspath(os.path.join(deepmd, os.pardir))
        elif dir_type == 'xyz':
            tmp = os.path.abspath(os.path.join(deepmd, os.pardir))
            recal_dir = os.path.abspath(os.path.join(tmp, os.pardir))
        outcar = os.path.join(recal_dir,'OUTCAR')
        fp = open(outcar)
        for line in fp:
            if 'nsw_tot' in line:
                print(line)
                nsw = int(line.split()[1])
                nsws.append(nsw)
                break
    return nsws
            
    
    
#fxyz = fxyz_recal_om8[1,2,4]
fxyz = [fxyz_recal_om8[0],fxyz_recal_om8[5]]
fxyz = [fxyz_recal_om8[0],fxyz_recal_pv[-1]]
fxyz = [fxyz_recal_om8[1],fxyz_recal_pv[0]]
fxyz = [fxyz_recal_om8[1],fxyz_recal_pv[0]]
#fxyz = fxyz_recal_pv[0]

fxyz = fxyz_recal_om8 + fxyz_recal_pv
nframes = get_nframes(fxyz,'xyz')


ele = 14 # Si

# make tags

if ele== 14 or ele == 12:
    n_atoms =  32
    n_frames = get_nframes(fxyz_recal_om8,'xyz')
    tag0 = sum(np.array(n_frames).astype(int))*n_atoms*[0]
    
    n_frames = get_nframes(fxyz_recal_pv,'xyz')
    tag1 = sum(np.array(n_frames).astype(int))*n_atoms*[1]
    tags = tag0 + tag1
#    tag0=(n_atoms*9 + n_atoms*250)*[0] + n_atoms*100*[1]
###########
########### The following is equivalent to above and have 
###########


asapxyz = ASAPXYZ(fxyz)


reduce_dict = {}
reduce_dict["preprocessing"] =  {"type": 'SCALE', 'parameter': None}
reduce_dict['skpca'] = {"type": 'SPARSE_KPCA',
                        'parameter':{"n_components": 3,
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

fig_spec = { 'outfile': 'test.png',
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
asap_plot.plot(proj[:, [0,1]],tags)



#asap_plot.plot(proj[:, [1,0]],plotcolor)

### plot command equivlaence
#plt.scatter(proj[:, [0]],proj[:, [1]],c = plotcolor)




