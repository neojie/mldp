#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 23:00:26 2020

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



#fxyz = fxyz_recal_om8[1,2,4]
fxyz = [fxyz_recal_om8[0],fxyz_recal_om8[5]]
fxyz = [fxyz_recal_om8[0],fxyz_recal_pv[-1]]
fxyz = [fxyz_recal_om8[1],fxyz_recal_pv[0]]
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

scale =True#False benchmark with the command line 
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
###########
###########
asapxyz = ASAPXYZ(fxyz)
reduce_dict = {}
reduce_dict['kpca'] = {"type": 'SPARSE_KPCA',
                        'parameter':{"n_components": 10,
                                     "n_sparse": -1, # no sparsification
                                "kernel": {"first_kernel": {"type": 'linear'}}}}

from asaplib.reducedim import Dimension_Reducers
dreducer = Dimension_Reducers(reduce_dict)


dm = asapxyz.fetch_computed_descriptors(['SOAP-n6-l6-c6.0-g0.44'])
asapxyz.fetch_computed_atomic_descriptors(['SOAP-n6-l6-c6.0-g0.44'])
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
plotcolor = np.ones(len(proj[:, [1,0]]))
asap_plot.plot(proj[:, [1,0]],plotcolor)
plotcolor = enthalpy[:]

asap_plot.plot(proj[:, [1,0]], plotcolor)







r1 = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r1-k111/recal/asap/map_per_atom.xyz'
r3 = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/map_per_atom.xyz'
fxyz = [r1,r3]

fxyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/merge_with_r1/map_per_atom.xyz'
fxyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/merge_with_r1/map_per_atom_ppv.xyz'

r1_tag = np.ones(100).astype(int)
r3_tag = (np.ones(250)*0).astype(int)
ppv_tag = (np.ones(75)*2).astype(int)
r1_r3_tag =np.concatenate((r1_tag,r3_tag))

r1_r3_ppv_tag =np.concatenate((r1_tag,r3_tag, ppv_tag))

np.savetxt('tag',r1_r3_tag,fmt='%d')
np.savetxt('tag_ppv',r1_r3_ppv_tag,fmt='%d')

#fmat = 'pca_coord'

fmat = 'skpca-d-10'
#fmat = '[*]'
#asapxyz = ASAPXYZ(fxyz)
a1 = ASAPXYZ(fxyz_recal_om8)
a2 = ASAPXYZ(fxyz_recal_pv)

asapxyz = ASAPXYZ(fxyz_recal_om8+fxyz_recal_pv)


dm, _ = asapxyz.get_descriptors(fmat, False)
dm_mg = asapxyz.get_atomic_descriptors(fmat, 12)
dm_oxygen = asapxyz.get_atomic_descriptors(fmat, 8)
dm_silicon = asapxyz.get_atomic_descriptors(fmat, 14)







plotcolor_volume, _, _, _ = set_color_function('volume', asapxyz)
plotcolor_density = np.zeros(len(plotcolor_volume))
for i in range(len(plotcolor_volume)):
    plotcolor_density[i] = 29.889703/plotcolor_volume[i]/3.
    
    
#tags = np.loadtxt('ice-54-labels.dat', dtype="str")[:,0]


#iceornot_hydrogen, _, _, _ = set_color_function('ice-or-not.tag', asapxyz, 0, 0, False, True, 1, False)
iceornot_oxygen, _, _, _ = set_color_function('tag', asapxyz, 0, 0, False, True, 8, False)
iceornot_silicon, _, _, _ = set_color_function('tag', asapxyz, 0, 0, False, True, 14, False)
iceornot_mg, _, _, _ = set_color_function('tag', asapxyz, 0, 0, False, True, 12, False)

iceornot_oxygen, _, _, _ = set_color_function('tag_ppv', asapxyz, 0, 0, False, True, 8, False)
iceornot_silicon, _, _, _ = set_color_function('tag_ppv', asapxyz, 0, 0, False, True, 14, False)
iceornot_mg, _, _, _ = set_color_function('tag_ppv', asapxyz, 0, 0, False, True, 12, False)

pc = [0, 1]
pca_d = len(pc)

# make plot
prefix = 'selected_phases'
fcolor = 'density'
plot_styles.set_nice_font()
fig, ax = plt.subplots()
fig.set_size_inches(10, 5)

cset1 = ax.scatter(dm[::-1, pc[0]],dm[::-1, pc[1]], c=plotcolor_density[::-1], 
           cmap='gnuplot', s=20, alpha=1.0, rasterized=True, label=None,
                                   vmax=None, vmin=None)

cbaxes = fig.add_axes([0.58, 0.85, 0.30, 0.02])
cbar=fig.colorbar(cset1, cax=cbaxes, orientation='horizontal')
cbar.ax.set_xlabel('Density [g/mL]',labelpad=-2,  size=14)
cbar.ax.tick_params(labelsize=14) 

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.tick_params(
    axis='x',       # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax.tick_params(
    axis='y',       # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off





plot_styles.set_nice_font()
fig, ax2 = plt.subplots()
fig.set_size_inches(10, 5)
ax2.scatter(dm_oxygen[::-1, 0], dm_oxygen[::-1, 1],c = list(range(len(dm_oxygen[::-1, 0]))),cmap='coolwarm')


plt.scatter(dm_silicon[::-1, 0], dm_silicon[::-1, 1],c = list(range(len(dm_silicon[::-1, 0]))),cmap='coolwarm')
plt.scatter(dm_mg[::-1, 0], dm_mg[::-1, 1], c = list(range(len(dm_mg[::-1, 0]))),cmap='coolwarm')



ax2.scatter(dm_oxygen[::-1, 0], dm_oxygen[::-1, 1], c=iceornot_oxygen[::-1], 
           cmap='coolwarm', s=2, alpha=1.0, rasterized=True,
                                   vmax=1, vmin=0)


plt.scatter(dm_oxygen[::-1, 0], dm_oxygen[::-1, 1], c=iceornot_oxygen[::-1], 
           cmap='winter', s=2, alpha=.3, rasterized=True,
                                   vmax=1, vmin=0)

plt.scatter(dm_silicon[::-1, 0], dm_silicon[::-1, 1], c=iceornot_silicon[::-1], 
           cmap='cividis', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.scatter(dm_mg[::-1, 0], dm_mg[::-1, 1], c=iceornot_mg[::-1], 
           cmap='winter', s=2, alpha=.03, rasterized=True,
                                   vmax=1, vmin=0)

plt.scatter(dm_mg[:, 0], dm_mg[:, 1], c=iceornot_mg, 
           cmap='coolwarm', s=2, alpha=.03, rasterized=True,
                                   vmax=1, vmin=0)

beg = 100*32
beg = 100*32
plt.scatter(dm_mg[:beg, 0], dm_mg[:beg, 1], c=iceornot_mg[:beg], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.xlim([5,23])
plt.ylim([-10,10])
plt.scatter(dm_mg[beg:, 0], dm_mg[beg:, 1], c=iceornot_mg[beg:], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)




plt.scatter(dm_silicon[beg:, 0], dm_silicon[beg:, 1], c=iceornot_silicon[beg:], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.xlim([10,28])
plt.ylim([-10,10])
plt.scatter(dm_silicon[:beg, 0], dm_silicon[:beg, 1], c=iceornot_silicon[:beg], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)


beg = 100*96

plt.scatter(dm_oxygen[beg:, 0], dm_oxygen[beg:, 1], c=iceornot_oxygen[beg:], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.xlim([-20,10])
plt.ylim([-25,30])
plt.scatter(dm_oxygen[:beg, 0], dm_oxygen[:beg, 1], c=iceornot_oxygen[:beg], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)






####

beg = 32*100
beg1 = 32*350
plt.xlim([5,30])
plt.ylim([-10,10])
plt.scatter(dm_silicon[:beg, 0], dm_silicon[:beg, 1], c=iceornot_silicon[:beg], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.scatter(dm_silicon[beg:beg1, 0], dm_silicon[beg:beg1, 1], c=iceornot_silicon[beg:beg1], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)
plt.xlim([5,30])
plt.ylim([-10,10])
plt.scatter(dm_silicon[beg1:, 0], dm_silicon[beg1:, 1], c=iceornot_silicon[beg1:], 
           cmap='coolwarm', s=2, alpha=.1, rasterized=True,
                                   vmax=1, vmin=0)

beg = 32*100
beg1 = 32*350
plt.xlim([5,20])
plt.ylim([-15,10])
plt.scatter(dm_mg[:beg, 0], dm_mg[:beg, 1], c=iceornot_mg[:beg], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)
plt.scatter(dm_mg[beg:beg1, 0], dm_mg[beg:beg1, 1], c=iceornot_mg[beg:beg1], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)
plt.xlim([5,20])
plt.ylim([-15,10])
plt.scatter(dm_mg[beg1:, 0], dm_mg[beg1:, 1], c=iceornot_mg[beg1:], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)



beg = 96*100
beg1 = 96*350
plt.xlim([-20,10])
plt.ylim([-20,20])
plt.scatter(dm_oxygen[:beg, 0], dm_oxygen[:beg, 1], c=iceornot_oxygen[:beg], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)
plt.scatter(dm_oxygen[beg:beg1, 0], dm_oxygen[beg:beg1, 1], c=iceornot_oxygen[beg:beg1], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)
plt.xlim([-20,10])
plt.ylim([-20,20])

plt.scatter(dm_oxygen[beg1:, 0], dm_oxygen[beg1:, 1], c=iceornot_oxygen[beg1:], 
           cmap='coolwarm', s=2, alpha=.2, rasterized=True,
                                   vmax=2, vmin=0)
