#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 18:48:56 2022

@author: jiedeng
"""
import matplotlib.pyplot as plt
import numpy as np


def render(inter,k):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    triangles = inter.triangulated_surface[0][inter.triangulated_surface[1]]
    ###
    pos   = inter.pos
    ### calculate distance
    verts, faces, normals = inter.triangulated_surface   
    # calculate proximity
    prox = k
    threshold = 0.7
    liq_in = prox>threshold
    vap_in = prox<threshold
        
    fig = plt.figure(figsize=(4, 5))
    ax1 = fig.add_subplot(111, projection='3d')
    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh1 = Poly3DCollection(triangles)
    mesh1.set_edgecolor('none')
    mesh1.set_alpha(0.3)
    # ax1.add_collection3d(mesh1)
    
    pos_liq = pos[liq_in]
    xyz_liq = np.vstack([pos_liq[::, 0], pos_liq[::, 1], pos_liq[::, 2]])
    
    ax1.scatter(xyz_liq[0],xyz_liq[1],xyz_liq[2],color='r')
    
    pos_vap = pos[vap_in]
    xyz_vap = np.vstack([pos_vap[::, 0], pos_vap[::, 1], pos_vap[::, 2]])
    
    ax1.scatter(xyz_vap[0],xyz_vap[1],xyz_vap[2],color='k')
    ax1.view_init(0, 0)
    
    ax1.set_xlabel("x-axis: a")
    ax1.set_ylabel("y-axis: b")
    ax1.set_zlabel("z-axis: c")
       
    plt.tight_layout()
    plt.show()

def plot_rho(z,rho,xlabel = r'$z (\AA)$', ylabel = r'$\rho (g/cm^{3})$' ):
    """
    plot density
    """
    fig,ax = plt.subplots(1,1)
    ax.plot(z,rho)
    ax.set_xlabel(xlabel,fontsize =13)
    ax.set_ylabel(ylabel,fontsize =13)
    fig.show()
    
    
    
def plot_ch(ch,name='stat.pdf'):
    # print
    fig,ax = plt.subplots(1,3,figsize=(16,4),sharey=False)
    ax[0].set_title('water content')
    ax[0].plot(ch[:,0],ch[:,1],label='solid')
    ax[0].plot(ch[:,0],ch[:,2],label='liquid')
    ax[0].plot(ch[:,0],ch[:,3],label='interface')
    ax[0].set_xlabel('step');ax[0].set_ylabel('water content')
    ax[0].legend()
    ax[1].set_title('scale of phases')
    ax[1].plot(ch[:,0],ch[:,4],label='solid')
    ax[1].plot(ch[:,0],ch[:,5],label='liquid')
    ax[1].plot(ch[:,0],ch[:,6],label='2w')
    ax[1].plot(ch[:,0],ch[:,6]/2,label='w')
    ax[1].legend()
    ax[1].set_xlabel('step');ax[1].set_ylabel('length (A)')
    ax[2].set_title('# of H')
    ax[2].plot(ch[:,0],ch[:,7],label='solid')
    ax[2].plot(ch[:,0],ch[:,8],label='liquid')
    ax[2].plot(ch[:,0],ch[:,9],label='2w')
    ax[2].legend()
    ax[2].set_xlabel('step');ax[2].set_ylabel('H count')
    fig.savefig(name,bbox_inches='tight')
