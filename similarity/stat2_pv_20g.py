#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 21:20:52 2022

@author: jiedeng
"""
import numpy as np
import ase.io
import MDAnalysis as mda
import matplotlib.pyplot as plt
from stat_lib import cal_inter, temporal_mass_density_proj_coarse, fit_gds_single_sided, fit_gds_double_sided


xyz     = '/Users/jiedeng/GD/papers/pv7_h/partition/data/pvw_20g_180pv9w_2600_cont1/merge.xyz'
project_axis = 2 #0,1,2 => x,y,z
solid_center = 11  ## do check!!! initial solid center, assuming solid center does not chagne drastically

#----------------------------------------------------------------------#

ase_xyz = ase.io.read(xyz,index=':') 
mda_xyz = mda.Universe(xyz)
length  =  len(ase_xyz)
num_atom = len(ase_xyz[0])
h_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==1]

length=1
ch = np.zeros((length,2))


def pbc_x(x,y,x0):
    """
    cut x-x0 portion of the x-y function and concatenate it in the end
    x[-1] is period
    """
    idx0 = np.argmin(np.abs(x-x0))
    ynew=np.concatenate((y[idx0:],y[0:idx0]))
    xnew=np.concatenate((x[idx0:],x[0:idx0]+x[-1]))
    return xnew,ynew

def unpbc_x(x0,dim):
    """
    undo pbc_x
    given x0, move it back to the original position given the dim dim
    """
    if x0 < dim:
        return x0
    else:
        return x0-dim

# def count_hs():
for fn in range(length):
    print('**'*50,fn)
    # build interface, mda_xyz does not contain cell dim info, ase_xyz cannot by used by pytim
    inter, k,ase_a,mda_a = cal_inter(ase_xyz,mda_xyz,fn,mesh=1., alpha=2.5, level=None)
    # build density map and project to target direction
    rho_zi,zi,rho_av     = temporal_mass_density_proj_coarse(inter,project_axis=project_axis,plot=False,level=None);
    zi_new,rho_zi_new    = pbc_x(zi,rho_zi,solid_center)
    dim                  = ase_a.get_cell_lengths_and_angles()[project_axis]
    # fit to double sided GDS function
    result,z0,z1,w = fit_gds_double_sided(zi_new,rho_zi_new,plot=True, verbose=False)
    # choose inteface to be w or 2w
    w=2*w
    # ensure z0 is close to solid center
    if (z0 - solid_center) <= (z1 - solid_center):
        pass
    else:
        tmp = z0
        z0 = z1
        z1 = tmp
    z1_unpbc = unpbc_x(z1,dim)
    
    ### solid z1, z1_unpbc
    """
    --------z0--------z1--
    unpbc
    --z1----z0------------
    
    or 
    ----z0-----z1---------
    -
    """
    # H positions in the projected axis
    hz = ase_a.positions[h_indices][:,project_axis]
    hz[hz<0]   += hz[hz<0] + dim
    hz[hz>dim] -= hz[hz>dim] - dim
        
    # H in solid positions in the projected axis 
    if z0 <  z1_unpbc:
        if z0 - w/2 >0:
            hs1 = hz[hz < (z0 - w/2)]
            zs1_bound = dim
        else:
            zs1_bound = z0 - w/2 + dim
            hs1 = np.array([])
        if z1_unpbc + w/2 < dim:
            tmp1 = hz[hz>(z1_unpbc + w/2)]
            hs2 = tmp1[tmp1<zs1_bound]
            zs2_bound = 0
        else:
            hs2 = np.array([])
            zs2_bound = z1_unpbc + w/2 - dim
        hs1   = hs1[hs1>zs2_bound]
        sol_z = np.concatenate((hs1,hs2))
    else:
        if (z1_unpbc+w/2) < (z0-w/2):
            tmp1=hz[hz < (z0 - w/2)]
            sol_z = tmp1[tmp1>(z1_unpbc + w/2)]
        else:
            sol_z = np.array([])

    if z0 <  z1_unpbc:
       if (z0+w/2) < (z1_unpbc-w/2):
           tmp1=hz[hz > (z0 + w/2)]
           liq_z = tmp1[tmp1<(z1_unpbc - w/2)]
       else:
           liq_z = np.array([])
    else:
        if z1_unpbc - w/2 >0:
            hl1 = hz[hz < (z1_unpbc - w/2)]
            zl1_bound = dim
        else:
            zl1_bound = z1_unpbc - w/2 + dim
            hl1 = np.array([])
        if z0 + w/2 < dim:
            tmp1 = hz[hz>(z0 + w/2)]
            hl2 = tmp1[tmp1<zl1_bound]
            zl2_bound = 0
        else:
            hl2 = np.array([])
            zl2_bound = z0 + w/2 -dim
        hl1 = hl1[hl1>zl2_bound]
        liq_z = np.concatenate((hl1,hl2))      
        
            
    
    # length of solid and liquid regime, interface thickness not considered
    if z0 <  z1_unpbc:
        ls,ly = dim -(z1_unpbc - z0) , (z1_unpbc - z0)
    else:
        ls, ly = z0-z1_unpbc, dim - (z0-z1_unpbc)
    print(sol_z,liq_z)
    # # of H in solid and liquid, respectively
    hs, hl = len(sol_z), len(liq_z)
    print(hs,hl)
    # concentration of h in solid and liquid, respectively
    chs, chl = hs/ls, hl/ly
    ch[fn,0],ch[fn,1] = chs, chl


plt.plot(zi,rho_zi,'o-')
plt.plot(zi_new,rho_zi_new,'*')



plt.figure()
plt.plot(ch[:fn,0])
plt.plot(ch[:fn,1])

# np.mean(ch[:fn,0])

np.mean(ch[:fn,0])/np.mean(ch[:fn,1])
np.argmax(ch[:,0])

np.nonzero(ch[:,0])
ch[np.nonzero(ch[:,0])]