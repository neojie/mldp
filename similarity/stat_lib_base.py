#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 10:47:24 2022

@author: jiedeng
"""
import numpy as np
# import ase.io
# import MDAnalysis as mda
import pytim
import matplotlib.pyplot as plt
from lmfit import Parameters, Model, report_fit
from scipy.interpolate import RegularGridInterpolator


params_global = None

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
    

def cal_inter(ase_xyz,mda_xyz,fn,mesh=1.5, alpha=2.5, level=None,mode='mass'):
    """
    calculate the interface with weight of atomic mass

    Parameters
    ----------
    ase_xyz : ase trajectories
    mda_xyz : mdanalysis trajectories
    fn : int, frame number
    mesh : float, optional
        DESCRIPTION. The default is 1.
    alpha : TYPE, optional
        DESCRIPTION. The default is 2.5.
    level : TYPE, optional
        DESCRIPTION. The default is None.
    mode:
        mass:  atomic mass as weighting factor
        k: k value as weighting factor

    Returns
    -------
    inter : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    ase_a : TYPE
        DESCRIPTION.
    mda_a : TYPE
        DESCRIPTION.

    """

    ase_a = ase_xyz[fn]
    mda_a = mda_xyz.trajectory[fn]
    mda_a.dimensions= ase_a.get_cell_lengths_and_angles()
    inter     = pytim.WillardChandler(mda_xyz, mesh=mesh, alpha=alpha,level=level)
    if mode == 'mass':
        k     = ase_a.get_masses()
    elif mode == 'k':
        k     = ase_a.arrays['k']
        ## re-calculate the density field using the fictious mass k
        inter.density_field,inter.mass_density_field = inter.kernel.evaluate_pbc_fast(inter.grid,k) 
        inter.mass_density_field_reshape = inter.mass_density_field.reshape(
            tuple(np.array(inter.ngrid).astype(int)))
    else:
        print('Weighting mode of {} is not supported!'.format(mode))
    _inter_density_two_end_equal(inter)
    inter.universe.atoms.masses = k   ## the intital masses for Si is problematic, bug in MDAnalysis
    return inter, k,ase_a,mda_a



    
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
    volume0 = volume/inter.ngrid[project_axis]  ## somehow treat Si as 0
    rho_av  = sum(inter.universe.atoms.masses)/volume    
    norm_factor = sum(inter.universe.atoms.masses)/volume0/proj.sum()
    return proj*norm_factor,rho_av


    
def temporal_mass_density_proj_coarse(inter,project_axis=1,plot=True,level=None):
    """
    plot coarse graining mass density profile file of ith step
    """
    rho_zi,rho_av     = mass_density_proj_coarse(inter, project_axis)
    zi = np.linspace(0,inter.box[project_axis],inter.ngrid[project_axis])
    if plot:
        from render import plot_1d_density
        plot_1d_density(zi,rho_zi)
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

def fit_gds_single_sided(zi,rho_zi,plot=True,verbose=True,vary_z0=True):
    """
    LS fitting to GDS, usually used for proximity vs. rho
    """
    params = Parameters()
    params.add('z0', value=0,vary=vary_z0)
    params.add('w' , value=1.5,vary=True)    
    params.add('rhov', value=min(rho_zi),min = 0.1,vary=True)
    params.add('rhol', value=max(rho_zi),min = 0,vary=True) 
    model = Model(gds_single_sided)
    result = model.fit(rho_zi,params,z = zi)
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
    dl = max(zi) - min(zi)
    z0, z1    = min(zi) + dl/6, max(zi) - dl/6
    rhol,rhov = min(rho_zi),max(rho_zi)
    zmin = zi[np.argmin(rho_zi)]
    
    params = Parameters()
    params.add('z0', value=z0,   min = min(zi),    max = zmin,   vary=True)
    params.add('z1', value=z1,   min = zmin,    max = max(zi),   vary=True)
    params.add('w' , value=3,    min = .1,     max =  max(zi),   vary=True)    
    params.add('rhov', value=rhov, min = rhol/4, max = rhov*4, vary=True)
    params.add('rhol', value=rhol, min = rhol/4, max = rhov*4, vary=True) 
    if params_global is None:
        pass
    else:
        params.add('z0', value=params_global.get('z0').value,   min = min(zi),    max = zmin,   vary=True)
        params.add('z1', value=params_global.get('z1').value,   min = zmin,    max = max(zi),   vary=True)
        params.add('w' , value=params_global.get('w').value,    min = .1,     max =  max(zi),   vary=True)    
        params.add('rhov', value=params_global.get('rhov').value, min = rhol/4, max = rhov*4, vary=True)
        params.add('rhol', value=params_global.get('rhol').value, min = rhol/4, max = rhov*4, vary=True)        
    model = Model(gds_double_sided,independent_vars=['z'])
    result = model.fit(rho_filter,params,z = zi_filter,weights=weights)

    if verbose:
        report_fit(result)
        result.params.pretty_print
    if plot:
        plt.figure()
        result.plot_fit()
    return result, result.params['z0'].value, result.params['z1'].value, result.params['w'].value


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

def count_z(z0,z1_unpbc, hz,w,dim,):
    """
    Given z coordinate, determine z that in liq and sol

    Parameters
    ----------
    z0 : TYPE
        DESCRIPTION.
    z1_unpbc : TYPE
        DESCRIPTION.
    hz : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    dim : TYPE
        DESCRIPTION.

    Returns
    -------
    sol_z : array
        DESCRIPTION.
    liq_z : array
        DESCRIPTION.

    """
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
    return sol_z, liq_z


    
def select_chi(ch,chi=0.2):
    ch1 = ch[ch[:,-1]<chi]
    return ch1

# def select_nonzero_sl(ch):
#     ch = ch[ch[:,4]>0.1]
#     ch = ch[ch[:,5]>0.1]
#     return ch

# def select_cs_less_cl(ch):
#     return ch[ch[:,1]<ch[:,2]]
    

    
############## proximity analysis


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
    

def plot_prox_rho(selected_prox, selected_rho,mean_prox, mean_rho,fitted_result):
    plt.plot(selected_prox, selected_rho,'.',alpha=0.1)
    plt.plot(mean_prox,mean_rho,'or')
    plt.plot(mean_prox,fitted_result.best_fit,'-r')
    plt.show()

    