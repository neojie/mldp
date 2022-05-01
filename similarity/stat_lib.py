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
    """
    calculate the interface with weight of k

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
    k     = ase_a.arrays['k']
    mda_a.dimensions= ase_a.get_cell_lengths_and_angles()
    inter     = pytim.WillardChandler(mda_xyz, mesh=mesh, alpha=alpha,level=level)
    inter.density_field,inter.mass_density_field = inter.kernel.evaluate_pbc_fast(inter.grid,k) # use similiarity as mass
    inter.mass_density_field_reshape = inter.mass_density_field.reshape(
        tuple(np.array(inter.ngrid).astype(int)))
    _inter_density_two_end_equal(inter)
    return inter, k,ase_a,mda_a

def cal_inter_mass(ase_xyz,mda_xyz,fn,mesh=1, alpha=2.5, level=None):
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
    k     = ase_a.get_masses()
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
    sol_z : TYPE
        DESCRIPTION.
    liq_z : TYPE
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

def analyze_single(idx,ase_xyz,mda_xyz, project_axis,solid_center0=-1,alpha = 2.5,plot=True,verbose=False, debug=False, assert_chi = True):
    """
    analyze all species

    Parameters
    ----------
    idx : TYPE
        DESCRIPTION.
    ase_xyz : TYPE
        DESCRIPTION.
    mda_xyz : TYPE
        DESCRIPTION.
    project_axis : TYPE
        DESCRIPTION.
    solid_center0 : TYPE, optional
        DESCRIPTION. The default is -1.
    alpha : TYPE, optional
        DESCRIPTION. The default is 2.5.
    plot : TYPE, optional
        DESCRIPTION. The default is True.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    debug : TYPE, optional
        DESCRIPTION. The default is False.
    assert_chi : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """
    print('**'*25,idx,'**'*25)
    num_atom = len(ase_xyz[0])
    h_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==1]
    mg_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==12]
    si_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==14]
    o_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==8]

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
    dim      = ase_a.get_cell_lengths_and_angles()[project_axis]
    lx,ly,lz = ase_a.get_cell_lengths_and_angles()[:3]
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

    mgz = ase_a.positions[mg_indices][:,project_axis]
    mgz[mgz<0]   += mgz[mgz<0] + dim
    mgz[mgz>dim] -= mgz[mgz>dim] - dim

    siz = ase_a.positions[si_indices][:,project_axis]
    siz[siz<0]   += siz[siz<0] + dim
    siz[siz>dim] -= siz[siz>dim] - dim

    oz = ase_a.positions[o_indices][:,project_axis]
    oz[oz<0]   += oz[oz<0] + dim
    oz[oz>dim] -= oz[oz>dim] - dim 
    
    sol_hz, liq_hz =count_z(z0,z1_unpbc, hz,w,dim)
    sol_mgz, liq_mgz =count_z(z0,z1_unpbc, mgz,w,dim)
    sol_siz, liq_siz =count_z(z0,z1_unpbc, siz,w,dim)
    sol_oz, liq_oz =count_z(z0,z1_unpbc, oz,w,dim)
    
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
    hs, hl = len(sol_hz), len(liq_hz)
    mgs, mgl = len(sol_mgz), len(liq_mgz)
    sis, sil = len(sol_siz), len(liq_siz)
    os, ol = len(sol_oz), len(liq_oz)

    # print(hs,hl)
    # concentration of h in solid and liquid, respectively
    hw = len(h_indices) - hs - hl
    mgw = len(mg_indices) - mgs - mgl
    siw = len(si_indices) - sis - sil
    ow = len(o_indices) - os - ol

    chs, chl, chw = hs/ls, hl/ll, (len(h_indices) - hs - hl)/2/w# interface width is 2w
    cmgs, cmgl, cmgw = mgs/ls, mgl/ll, (len(mg_indices) - mgs - mgl)/2/w# interface width is 2w
    csis, csil, csiw = sis/ls, sil/ll, (len(si_indices) - sis - sil)/2/w
    cos, col, cow = os/ls, ol/ll, (len(o_indices) - os - ol)/2/w    
    return chs, chl,chw, ls, ll, w,hs, hl , hw, z0, z1_unpbc, result.redchi, lx,ly,lz, mgs,mgl,mgw,sis,sil,siw,os,ol,ow    
    
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
    ch = np.zeros((size,25))  # 
    for fn in range(start_idx,en_idx,step):
        chs, chl,chw, ls, ll, w, hs, hl, hw,z0,z1_unpbc,chi,lx,ly,lz,mgs,mgl,mgw,sis,sil,siw,os,ol,ow = analyze_single(fn,ase_xyz,mda_xyz, project_axis,alpha=alpha,plot=False,assert_chi=assert_chi)
        idx = (fn - start_idx)//step
        ch[idx,0] = fn
        ch[idx,1], ch[idx,2], ch[idx,3] = chs, chl, chw
        ch[idx,4], ch[idx,5], ch[idx,6] = ls, ll, w
        ch[idx,7], ch[idx,8], ch[idx,9] = hs, hl , hw
        ch[idx,10], ch[idx,11] = z0,z1_unpbc
        ch[idx,12] = chi
        ch[idx,13],ch[idx,14], ch[idx,15] = lx,ly,lz 
        ch[idx,16],ch[idx,17], ch[idx,18] = mgs,mgl,mgw 
        ch[idx,19],ch[idx,20], ch[idx,21] = sis,sil,siw
        ch[idx,22],ch[idx,23], ch[idx,24] = os,ol,ow

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
    fmt = '%d %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d %d %d %d %d %d %d'
    np.savetxt(name,ch,fmt=fmt,
                header='{2} \n python start_idx = {0} alpha = {1} dim = ls+ll+w*2 \n id,cs,cl,cw,ls,ll.lw,hs,hl,hw,z0,z1_unpbc,reduced_chi, lx,ly,lz, mgs,mgl,mgw,sis,sil,siw,os,ol,ow'.format(start_idx, alpha, filepath))   

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


    
def select_chi(ch,chi=0.2):
    ch1 = ch[ch[:,12]<chi]
    return ch1

def select_nonzero_sl(ch):
    ch = ch[ch[:,4]>0.1]
    ch = ch[ch[:,5]>0.1]
    return ch

def select_cs_less_cl(ch):
    return ch[ch[:,1]<ch[:,2]]

def show_atoms_num(ch):
    fig,ax = plt.subplots(1,3,figsize=(15,4),sharey=True)
    ax[0].set_title('solid')
    ax[0].plot(ch[:,0],ch[:,16],label='Mg = {0:.1f}'.format(ch[:,16].mean()))
    ax[0].plot(ch[:,0],ch[:,19],label='Si = {0:.1f}'.format(ch[:,19].mean()))
    ax[0].plot(ch[:,0],ch[:,22],label='O = {0:.1f}'.format(ch[:,22].mean()))
    ax[0].plot(ch[:,0],ch[:,7],label='H = {0:.1f}'.format(ch[:,7].mean()))
    ax[0].set_xlabel('step');ax[0].set_ylabel('# of atoms')
    ax[0].legend()
    ax[1].set_title('liquid')
    ax[1].plot(ch[:,0],ch[:,17],label='Mg = {0:.1f}'.format(ch[:,17].mean()))
    ax[1].plot(ch[:,0],ch[:,20],label='Si =  {0:.1f}'.format(ch[:,20].mean()))
    ax[1].plot(ch[:,0],ch[:,23],label='O = {0:.1f}'.format(ch[:,23].mean()))
    ax[1].plot(ch[:,0],ch[:,8],label='H = {0:.1f}'.format(ch[:,8].mean()))
    ax[1].legend()
    ax[1].set_xlabel('step');
    ax[2].set_title('interface')
    ax[2].plot(ch[:,0],ch[:,18],label='Mg= {0:.1f}'.format(ch[:,18].mean()))
    ax[2].plot(ch[:,0],ch[:,21],label='Si= {0:.1f}'.format(ch[:,21].mean()))
    ax[2].plot(ch[:,0],ch[:,24],label='O= {0:.1f}'.format(ch[:,24].mean()))
    ax[2].plot(ch[:,0],ch[:,9],label='H = {0:.1f}'.format(ch[:,9].mean()))
    ax[2].legend()
    ax[2].set_xlabel('step');
    fig.savefig('atoms.pdf',bbox_inches='tight')
    
def show_water_step(ch,water):
    w2step = []
    for i in range(len(ch)):
        w2step.append(np.mean(water[:i]))
    plt.figure()
    plt.plot(ch[:,0],w2step)
    plt.xlabel('step')
    plt.ylabel('D=Cs/Cl (Xh2O/(Xh2O+XSi))')
    plt.grid()
    plt.savefig('stepvsD.pdf',bbox_inches='tight')
    return w2step

def show_length_scale(ch):
    fig,ax = plt.subplots(1,3,figsize=(15,4),sharey=False)
    ax[0].plot(ch[:,0],ch[:,4],label='s = {0:.1f}'.format(ch[:,4].mean()))
    ax[0].plot(ch[:,0],ch[:,5],label='l = {0:.1f}'.format(ch[:,5].mean()))
    ax[0].plot(ch[:,0],ch[:,6]*2,label='4w = {0:.1f}'.format(ch[:,6].mean()*2))
    ax[0].set_xlabel('step');
    ax[0].legend()
    ax[0].set_title('length (A)')
    ax[1].set_title('volume (A^3)')
    vs = ch[:,4]*ch[:,13]*ch[:,14]
    vl = ch[:,5]*ch[:,13]*ch[:,14]
    vw = ch[:,6]*2*ch[:,13]*ch[:,14]
    ax[1].plot(ch[:,0],vs,label='s = {0:.1f}'.format(vs.mean()))
    ax[1].plot(ch[:,0],vl,label='l =  {0:.1f}'.format(vl.mean()))
    ax[1].plot(ch[:,0],vw,label='4w = {0:.1f}'.format(vw.mean()))
    ax[1].legend()
    ax[1].set_xlabel('step'); ax[1].set_ylabel('volume (A^3)')
    
    gmolA2gcm = 1/6.0221409e23/(1e-30)/1e6#
    mmg, mgsi, mgo, mgh = 24.305, 28.085, 15.999, 1.008
    rhos=(ch[:,16]*mmg+ch[:,19]*mgsi+ch[:,22]*mgo+ch[:,7]*mgh)/vs*gmolA2gcm
    rhol=(ch[:,17]*mmg+ch[:,20]*mgsi+ch[:,23]*mgo+ch[:,8]*mgh)/vl*gmolA2gcm
    rhow = (ch[:,18]*mmg+ch[:,21]*mgsi+ch[:,24]*mgo+ch[:,9]*mgh)/vw*gmolA2gcm

    ax[2].plot(ch[:,0],rhos,label='s = {0:.1f}'.format(rhos.mean()))
    ax[2].plot(ch[:,0],rhol,label='l =  {0:.1f}'.format(rhol.mean()))
    ax[2].plot(ch[:,0],rhow,alpha=0.5,label='4w = {0:.1f}'.format(rhow.mean()))
    ax[2].legend()
    ax[2].set_xlabel('step'); ax[2].set_title('rho (g/cm^3)')
    fig.savefig('scale.pdf',bbox_inches='tight')
    return vs, vl, vw, rhos, rhol, rhow

def show_water_content(ch):
    xs=ch[:,7]/2/(ch[:,19]+ch[:,7]/2)  # water molal fraction mol of water/(mole of H/2 + mol of Si)
    xl=ch[:,8]/2/(ch[:,20]+ch[:,8]/2)
    xw=ch[:,9]/2/(ch[:,21]+ch[:,9]/2)
    fig,ax = plt.subplots(1,1,figsize=(5,4))
    ax.plot(ch[:,0],xs,label='s = {0:.5f}'.format(xs.mean()))
    ax.plot(ch[:,0],xl,label='l = {0:.5f}'.format(xl.mean()))
    ax.plot(ch[:,0],xw,label='4w = {0:.5f}'.format(xw.mean()))
    ax.set_xlabel('step');ax.set_ylabel('water content' )
    ax.legend()
    fig.savefig('water.pdf',bbox_inches='tight')
    return xs, xl, xw
       
def show(ch,beg=0,end=-1):
    ch = ch[beg:end]
    print('++'*20)
    print(' '*10 + 'Do select frames based on chi, and nonzero_sl' + ' '*10)
    print('++'*20)
    print(' '*10 + 'Step 1: # of atoms  in each phase vs. step' + ' '*10)
    show_atoms_num(ch)
    
    print(' '*10 + 'Step 2: z length, volume, and density of each phase' + ' '*10)
    vs, vl, vw, rhos, rhol, rhow = show_length_scale(ch)
    
    print(' '*10 + 'Step 3: water content in each phase (H/2)/(Si+H/2)' + ' '*10)
    xs, xl, xw = show_water_content(ch)
    
    print(' '*10 + 'Step 4: water partitioning vs. step' + ' '*10)
    w2step = show_water_step(ch,xs/xl)
    print('++'*20)
    
    print(' '*10 + 'Step 5: Print out statistics' + ' '*10)
    
    print('solid: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,16].mean()*2 + ch[:,19].mean()*4 + ch[:,7].mean(), 2*ch[:,22].mean()))
    print('liquid: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,17].mean()*2 + ch[:,20].mean()*4 + ch[:,8].mean(), 2*ch[:,23].mean()))
    print('interface: 2Mg+4Si+1H = {0}, 2O = {1}'.format(ch[:,18].mean()*2 + ch[:,21].mean()*4 + ch[:,9].mean(), 2*ch[:,24].mean()))

    print("solid (A^3):", vs.mean())
    print("liquid (A^3):", vl.mean())
    print("interface (A^3):", vw.mean())

    print("solid xs:, water content(ppm)", xs.mean(), xs.mean()*18/(xs.mean()*18+100)*1e6 ) # add block average with errorbar
    print("liquid xs:, liquid content(ppm)", xl.mean(), xl.mean()*18/(xl.mean()*18+100)*1e6) # add block average with errorbar
    print("D", xs.mean()/xl.mean(), w2step[-1])
    