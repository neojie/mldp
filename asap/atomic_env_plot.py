#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 23:00:26 2020

@author: jiedeng
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from asaplib.data import ASAPXYZ
#from asaplib.reducedim import PCA
from asaplib.plot import *

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
asapxyz = ASAPXYZ(fxyz)

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
