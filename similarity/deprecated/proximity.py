#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 14:09:55 2022

@author: jiedeng
"""
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
from stat_lib_base import fit_gds_single_sided

def cal_proximity(inter):
    """
    calculate proximity
    """
    pos   = inter.original_positions
    ### calculate distance    
    verts, faces, normals = inter.triangulated_surface  
    proximity = np.zeros(len(pos))
    tag = np.zeros(len(pos))  
    for num, point in enumerate(pos):
        RS = (point - verts) # R to surface vector
        tmp = RS**2
        dis = np.sum(tmp,axis=1)**.5
        min_ind = np.argmin(np.abs(dis))      # min distance index for the vertex
        min_nor = normals[min_ind]            # normal vector corresponding to min distance
        min_dis = dis[min_ind]
        tmp2    = np.sum(min_nor*RS[min_ind])
        min_dis = -tmp2/np.abs(tmp2)*min_dis # same direction, out of the interface
        proximity[num] = min_dis              # based on test
        if normals[min_ind,2]>0:
            tag[num] = 1                #upper interface
        else:
            tag[num] = -1
    return proximity,tag


def cal_local_rho(inter):
    x = np.unique(inter.grid[0])
    y = np.unique(inter.grid[1])
    z = np.unique(inter.grid[2])
    V = inter.mass_density_field_reshape
    # mass_density_field_reshape shape is ngrid[0]*ngrid[1]*ngrid[2]
    fn = RegularGridInterpolator([x,y,z], V)
    out_rho = []
    selected_indice = []
    # out_prox = []
    for num,pos in enumerate(inter.original_positions):
        try:
            out_rho.append(fn(pos))
            selected_indice.append(num)
            # out_prox.append(proximity[num])
        except:
            pass
    out_rho =  np.array(out_rho)
    return out_rho, selected_indice

    
def plot_3d_interface(inter):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    triangles = inter.triangulated_surface[0][inter.triangulated_surface[1]]
    ###
    pos   = inter.original_positions
    ### calculate distance
    verts, faces, normals = inter.triangulated_surface   
    # calculate proximity
    prox, tag = cal_proximity(inter)
    
    liq_in = prox>0
    vap_in = prox<0
    
    fig = plt.figure(figsize=(4, 5))
    ax1 = fig.add_subplot(111, projection='3d')
    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh1 = Poly3DCollection(triangles)
    mesh1.set_edgecolor('none')
    mesh1.set_alpha(0.3)
    ax1.add_collection3d(mesh1)
    
    pos_liq = pos[liq_in]
    xyz_liq = np.vstack([pos_liq[::, 0], pos_liq[::, 1], pos_liq[::, 2]])
    ax1.scatter(xyz_liq[0],xyz_liq[1],xyz_liq[2],color='r')
    
    pos_vap = pos[vap_in]
    xyz_vap = np.vstack([pos_vap[::, 0], pos_vap[::, 1], pos_vap[::, 2]])
    
    ax1.scatter(xyz_vap[0],xyz_vap[1],xyz_vap[2],color='w')
    # ax1.view_init(0, 135)
    for ii in range(0,360,40):
        ax1.view_init(elev=10., azim=ii)
        # plt.savefig("movie%d.png" % ii)   
    ax1.set_xlabel("x-axis: a")
    ax1.set_ylabel("y-axis: b")
    ax1.set_zlabel("z-axis: c")
   
    plt.tight_layout()
    plt.show()
    

def _average_prox_vs_rho(prox,rho,nbins):
    prox_new = np.linspace(min(prox),max(prox),nbins)
    rho_new = np.zeros(nbins-1)
    for ii in range(nbins-1):
        rho_new[ii] = np.mean(rho[(prox>prox_new[ii]) & (prox<prox_new[ii+1])])
    return prox_new[:-1],rho_new
    

def plot_prox_rho(inter,proximity,fitted_result):

    plt.plot(selected_prox, selected_rho,'.',alpha=0.1)
    plt.plot(mean_prox,mean_rho,'or')
    plt.plot(mean_prox,fitted_result.best_fit,'-r')
    plt.show()
    # average over proximity
    # coord=[]
    # for i in range(3):
    #     coord.append(np.linspace(0,inter.box[i],inter.ngrid[i]))
    # fn = RegularGridInterpolator(coord, V)
    # pts = pos
    # print(fn(pts))


