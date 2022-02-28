#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 10:47:24 2022

@author: jiedeng
"""
import numpy as np
import ase.io
import MDAnalysis as mda
import pytim
import matplotlib.pyplot as plt
from lmfit import Parameters, Model, report_fit

def _inter_density_two_end_equal(inter):
    """
    the 0, and end density does not match, probably due to some corner case
    if a dege is cut evenly into n pieces, there is n+1 grid point,
    if we only calculate the density in the center of the n pieces then the end
    and start would not match. But this does not seem to the case
    
    /Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/pytim/gaussian_kde_pbc.py
    
    """
    
    rho = inter.mass_density_field_reshape

    ave = (rho[0,:,:]+rho[-1,:,:])/2
    rho[0,:,:] = rho[-1,:,:] = ave
    
    ave = (rho[:,0,:]+rho[:,-1,:])/2
    rho[:,0,:] = rho[:,-1,:] = ave
    
    ave = (rho[:,:,0]+rho[:,:,-1])/2
    rho[:,:,0] = rho[:,:,-1] = ave
    inter.mass_density_field_reshape = rho
    
def cal_inter(ase_xyz,mda_xyz,fn,mesh=1, alpha=2.5, level=None):
    ase_a = ase_xyz[fn]
    mda_a = mda_xyz.trajectory[fn]
    k     = ase_a.arrays['k']
    mda_a.dimensions= ase_a.get_cell_lengths_and_angles()
    inter     = pytim.WillardChandler(mda_xyz, mesh=mesh, alpha=alpha,level=level)
    inter.density_field,inter.mass_density_field = inter.kernel.evaluate_pbc_fast(inter.grid,k) # use similiarity as mass
    inter.mass_density_field_reshape = inter.mass_density_field.reshape(
        tuple(np.array(inter.ngrid).astype(int)))
    _inter_density_two_end_equal(inter)
    return inter, k,ase_a,mda_a



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
    
# density_field,mass_density_field = inter.kernel.evaluate_pbc_fast(inter.grid,inter.universe.atoms.masses)


def projection(matrix, project_axis):
    """
    project a 3D matrix to x,y,z direction
    matrix : a 3D numpy array
    project_axis : int, 1,2,3  => x,y,z
    """
    rules = {0:(1,1),1:(0,1),2:(0,0)}
    i,j=rules[project_axis]
    return np.sum(np.sum(matrix,axis=i),axis=j)
    
def mass_density_proj_coarse(inter,project_axis):
    """
    Calculate mass density profile along z with coarse graining
    """
    NA      = 6.022e23
    dens_re = inter.mass_density_field_reshape
    proj    = projection(dens_re,project_axis)
    volume  = inter.box[0]*inter.box[1]*inter.box[2]*NA*1e-24
    volume0 = volume/inter.ngrid[-1]
    rho_av  = sum(inter.universe.atoms.masses)/volume    
    norm_factor = sum(inter.universe.atoms.masses)/volume0/proj.sum()
    return proj*norm_factor,rho_av

def plot_rho(z,rho,xlabel = r'$z (\AA)$', ylabel = r'$\rho (g/cm^{3})$' ):
    """
    plot density
    """
    fig,ax = plt.subplots(1,1)
    ax.plot(z,rho)
    ax.set_xlabel(xlabel,fontsize =13)
    ax.set_ylabel(ylabel,fontsize =13)
    fig.show()
    
def temporal_mass_density_proj_coarse(inter,project_axis=1,plot=True,level=None):
    """
    plot coarse graining mass density profile file of ith step
    """
    rho_zi,rho_av     = mass_density_proj_coarse(inter, project_axis)
    zi = np.linspace(0,inter.box[project_axis],inter.ngrid[project_axis])
    if plot:
        plot_rho(zi,rho_zi)
    return rho_zi,zi,rho_av


def gds_single_sided(z,w,rhol,rhov,z0=0):
    """
    Gibbs dividing surface
    ---------
    input:
        rhol: liquid
        rhov: vapor
        w: width of the dividing surface
    output:
        rho:    
    """
    rho = .5*(rhol + rhov) + .5*(rhol-rhov)*np.tanh((z-z0)/w)
    return rho

def fit_gds_single_sided(zi,rho_zi,plot=True,verbose=True):
    """
    LS fitting to GDS
    """
    zi_filter  = zi[rho_zi>0]
    rho_filter = rho_zi[rho_zi>0]
    params = Parameters()
    params.add('z0', value=15,vary=True)
    params.add('w' , value=1.5,vary=True)    
    params.add('rhov', value=0.2,min = 0.1,vary=True)
    params.add('rhol', value=2.3,min = 0,vary=True) 
    model = Model(gds_single_sided)
    half = int(len(rho_filter))
    result = model.fit(rho_filter[:half],params,z = zi_filter[:half])
    if verbose:
        report_fit(result)
        result.params.pretty_print
    if plot:
        result.plot_fit()
    return result, result.params['z0'].value, result.params['w'].value


def gds_double_sided(z,w,rhol,rhov,z0=0,z1=5):
    """
    Gibbs dividing surface
    ---------
    input:
        rhol: liquid
        rhov: vapor
        w: width of the dividing surface
    output:
        rho:    
    """
#    rho = .5*(rhol + rhov) + .5*(rhol-rhov)*(np.tanh((z-z0)/w) - np.tanh((z-z1)/w))
    rho = rhov + .5*(rhol-rhov)*(np.tanh((z-z0)/w) - np.tanh((z-z1)/w))

    return rho

def fit_gds_double_sided(zi,rho_zi,weights=None,
                         z0min = 0, z1min = 0,
                         z0max = 400, z1max = 400,
                         rhovmin = 0, rholmin = 0,
                         rhovmax = 10, rholmax = 10,
                         wmin = 0, wmax = 100,
                         plot=True, verbose=True):
    """
    LS fitting to GDS
    """
#    zi_filter  = zi[rho_zi>0]
#    rho_filter = rho_zi[rho_zi>0]
    zi_filter  = zi
    rho_filter = rho_zi
    z0, z1 = max(zi)/4, max(zi)*3/4
    rhol,rhov = min(rho_zi),max(rho_zi)
    params = Parameters()
    params.add('z0', value=z0,   min = min(zi),    max = max(zi),   vary=True)
    params.add('z1', value=z1,   min = min(zi),    max = max(zi),   vary=True)
    params.add('w' , value=5,    min = 1,     max =  max(zi)/2,   vary=True)    
    params.add('rhov', value=rhol, min = rhovmin, max = rhovmax, vary=True)
    params.add('rhol', value=rhov, min = rholmin, max = rholmax, vary=True) 
    model = Model(gds_double_sided)
    result = model.fit(rho_filter,params,z = zi_filter,weights=weights)
    assert result.redchi<0.2
    if verbose:
        report_fit(result)
        result.params.pretty_print
    if plot:
        result.plot_fit()
    return result, result.params['z0'].value, result.params['z1'].value, result.params['w'].value
# inter.kernel.evaluate_pbc_fast()
# plot_temporal_inter()