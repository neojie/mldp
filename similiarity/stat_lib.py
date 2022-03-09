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
                         plot=True, verbose=True,debug=False):
    """
    LS fitting to GDS
    """
#    zi_filter  = zi[rho_zi>0]
#    rho_filter = rho_zi[rho_zi>0]
    zi_filter  = zi
    rho_filter = rho_zi
    z0, z1    = max(zi)/4, max(zi)*3/4
    rhol,rhov = min(rho_zi),max(rho_zi)
    zmin = zi[np.argmin(rho_zi)]
    
    params = Parameters()
    params.add('z0', value=z0,   min = min(zi),    max = zmin,   vary=True)
    params.add('z1', value=z1,   min = zmin,    max = max(zi),   vary=True)
    params.add('w' , value=5,    min = 1,     max =  max(zi)/2,   vary=True)    
    params.add('rhov', value=rhov, min = rhol/1.2, max = rhov*1.2, vary=True)
    params.add('rhol', value=rhol, min = rhol/1.2, max = rhov*1.2, vary=True) 
    model = Model(gds_double_sided)
    result = model.fit(rho_filter,params,z = zi_filter,weights=weights)

    if verbose:
        report_fit(result)
        result.params.pretty_print
    if plot:
        plt.figure()
        result.plot_fit()
    return result, result.params['z0'].value, result.params['z1'].value, result.params['w'].value
# inter.kernel.evaluate_pbc_fast()
# plot_temporal_inter()


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

def analyze_single(idx,ase_xyz,mda_xyz, project_axis,solid_center0=-1,alpha = 2.5,plot=True,verbose=False, debug=False, assert_chi = True):
    print('**'*25,idx,'**'*25)
    num_atom = len(ase_xyz[0])
    h_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==1]
    # build interface, mda_xyz does not contain cell dim info, ase_xyz cannot by used by pytim
    inter, k,ase_a,mda_a = cal_inter(ase_xyz,mda_xyz,idx,mesh=1., alpha=alpha, level=None)
    # build density map and project to target direction
    rho_zi,zi,rho_av     = temporal_mass_density_proj_coarse(inter,project_axis=project_axis,plot=plot,level=None);
    # find solid_center
    if solid_center0 <0: # if solid center is not specified, use the default value
        solid_center0 = zi[np.argmax(rho_zi)]
    for  solid_center in [solid_center0, solid_center0 -2, solid_center0+2, solid_center0-1, solid_center0+1]:
        zi_new,rho_zi_new    = pbc_x(zi,rho_zi,solid_center)
        # fit to double sided GDS function
        result,z0,z1,w1 = fit_gds_double_sided(zi_new,rho_zi_new,plot=plot, verbose=verbose,debug=debug)
        if result.redchi<0.2:
            break
    if debug:
        print("solid center, initial solid center ::", solid_center, solid_center0)
    if assert_chi:
        assert result.redchi<0.5        
    # choose inteface to be w or 2w
    w=2*w1
    # ensure z0 is close to solid center
    if (z0 - solid_center) <= (z1 - solid_center):
        pass
    else:
        tmp = z0
        z0 = z1
        z1 = tmp
    dim                  = ase_a.get_cell_lengths_and_angles()[project_axis]
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
        ls,ll = dim -(z1_unpbc - z0) - w, (z1_unpbc - z0) - w
    else:
        ls, ll = z0-z1_unpbc - w, dim - (z0-z1_unpbc) - w
    
    if ls <0:
        ls = 0.00001
    if ll <0:
        ll = 0.00001
    # print(sol_z,liq_z)
    # # of H in solid and liquid, respectively
    hs, hl = len(sol_z), len(liq_z)
    # print(hs,hl)
    # concentration of h in solid and liquid, respectively
    hw = len(h_indices) - hs - hl
    chs, chl, chw = hs/ls, hl/ll, (len(h_indices) - hs - hl)/2/w# interface width is 2w
    
    return chs, chl,chw, ls, ll, w,hs, hl , hw, z0, z1_unpbc, result.redchi

def analyze(start_idx,en_idx, filepath,ase_xyz,mda_xyz, project_axis,alpha, step=1,save=True,name='stat.txt', assert_chi=True):
    """
    analyze multiple frames 

    Parameters
    ----------
    start_index : int, 
        start index, count from 0.
    length : int
        num of frames to be analyzed.
    save :  bool, optional
        save output,. The default is True.
    name : string, optional
        output name. The default is 'stat.txt'.

    Returns
    -------
    None.

    """
    size = len(range(start_idx,en_idx,step))
    ch = np.zeros((size,13))  # 
    for fn in range(start_idx,en_idx,step):
        chs, chl,chw, ls, ll, w, hs, hl, hw,z0,z1_unpbc,chi = analyze_single(fn,ase_xyz,mda_xyz, project_axis,alpha=alpha,plot=False,assert_chi=assert_chi)
        idx = (fn - start_idx)//step
        ch[idx,0] = fn
        ch[idx,1], ch[idx,2], ch[idx,3] = chs, chl, chw
        ch[idx,4], ch[idx,5], ch[idx,6] = ls, ll, w
        ch[idx,7], ch[idx,8], ch[idx,9] = hs, hl , hw
        ch[idx,10], ch[idx,11] = z0,z1_unpbc
        ch[idx,12] = chi
    if save:
        save_ch(ch, start_idx, alpha, filepath,name)
    return ch

def save_ch(ch,start_idx, alpha, filepath,name):
    """
    # cs (H content in solid), cl (H content in liquid), cw (H in interface)
    # ls (solid slab length),ll (liquid slab length),w (interface), 
    # hs (# of H in solid), hl (# of H in liquid),hw (# of H in inteface)

    Parameters
    ----------
    ch : TYPE
        DESCRIPTION.
    start_idx : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    filepath : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fmt = '%d %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d %.4f %.4f %.4f'
    np.savetxt(name,ch,fmt=fmt,
                header='{2} \n python start_idx = {0} alpha = {1} dim = ls+ll+w*2 \n id,cs,cl,cw,ls,ll.lw,hs,hl,hw,z0,z1_unpbc,reduced_chi'.format(start_idx, alpha, filepath))    

def show(ch):
    # print
    fig,ax = plt.subplots(1,2,figsize=(12,4),sharey=False)
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
    fig.savefig('stat.pdf',bbox_inches='tight')
    result1=np.mean(ch[:,1])/np.mean(ch[:,2])
    result2=np.mean(ch[:,1][np.nonzero(ch[:,2])]/ch[:,2][np.nonzero(ch[:,2])])
    print("**"*25+'water in solid / water in liquid'+"**"*25)
    print("cs/cl",result1, result2, "(recommened)")
    print("**"*25+'frames with non-zero H in solid'+"**"*25)
    print(ch[:,0][np.nonzero(ch[:,1])])
    
