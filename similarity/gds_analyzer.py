#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 19:07:16 2022
TODO 
1- proximity for irregular plane
each frame should have one inter object
each inter has its own ele masses, etc, not analyze
analyzer is high level that only controls the analysis setup
@author: jiedeng
"""
import numpy as np
import ase.io
import MDAnalysis as mda
from stat_lib_base import cal_inter, fit_gds_double_sided,  temporal_mass_density_proj_coarse, pbc_x, unpbc_x, count_z, \
    cal_proximity, cal_local_rho, plot_3d_interface, plot_prox_rho, _average_prox_vs_rho, fit_gds_single_sided, cal_mol_frac, _unique_unsorted
import matplotlib.pyplot as plt

class GDSAnalyzer(object):
    def __init__(self,
             xyz  = 'merge.xyz',
             mode = 'mass',
             nw   = 2, # 
             project_axis = 2,
             begin = 0,
             end   = -1,
             step  = 1,
             alpha = 2.5,
             verbose  = False,
             is_analyze_planar = True,
             is_analyze_proximity = True,
             is_analyze_mol_frac = True):
        self.xyz = xyz
        self.ase_xyz = ase.io.read(xyz,index=':') 
        self.mda_xyz = mda.Universe(xyz)
        self.mode = mode
        self.nw = nw
        self.begin =begin
        if end == -1:
            self.end = len(self.mda_xyz.trajectory)
        else:
            self.end = end
        self.step = step
        self.project_axis = project_axis
        self.alpha = alpha
        self.mesh = 1
        self.level = None
        self.verbose  = verbose
        self.is_analyze_mol_frac  = is_analyze_mol_frac
        self.is_analyze_proximity = is_analyze_proximity
        self.is_analyze_planar    = is_analyze_planar

        if self.is_analyze_mol_frac:
            self.mol_frac_frame_count = 0
        
        self.num_atom = len(self.mda_xyz.trajectory[0])
        
        self.unique_ele  = _unique_unsorted(self.mda_xyz.universe.atoms.elements)
        self.unique_mass = _unique_unsorted(self.mda_xyz.universe.atoms.masses)
        
        self.num_ele  = len(self.unique_ele)
        
        self.ele_atomic_number = np.unique(self.ase_xyz[0].get_atomic_numbers())
        # ele_atomic_number = np.unique(tmp.get_atomic_numbers())

        self.ele_chemical_symbol = []
        for num in self.ele_atomic_number:
            self.ele_chemical_symbol.append(ase.data.chemical_symbols[num])
        # print('System:', self.ase_xyz[0].get_chemical_formula())
        
        self.ele_indices = []
        for num in self.ele_atomic_number:
            self.ele_indices.append(np.array(range(self.num_atom))[self.ase_xyz[0].get_atomic_numbers()==num])
        # h_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==2]  #2
        # mg_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==12] # mg
        # si_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==26] # fe
        # o_indices      = np.array(range(num_atom))[ase_xyz[0].get_atomic_numbers()==8]  # o
        self.sum_elements_counts = []
        self.sum_phase_dimension = []
        self.sum_rechi =[]
        self.sum_prox = []
        if self.is_analyze_planar or self.is_analyze_proximity:
            self.analyze()
            self.save_sum()            
        # if analyze_planar:
        #     self.analyze()
        #     self.save_sum()
        # if analyze_proximity:
        #     self.analyze_prox()
        
    def save_sum(self):
        if self.is_analyze_planar:
            name = 'sum_counts_{0}_{1}.txt'.format(self.begin,self.end-1)
            fmt =  '%d    ' + '%d %d %d    '*len(self.ele_atomic_number) +'%.4f' # for each element, we have three phases
            dat = np.concatenate((np.array([range(self.begin,self.end)]).T, self.sum_elements_counts,np.array([self.sum_rechi]).T),axis=1)
            np.savetxt(name,dat,fmt=fmt,
                        header='{0} \n {3} \n {1} \n {2}'.format(self.xyz,'solid liquid interface','id '+'           '.join(self.ele_chemical_symbol) +'    chi', 'nw = {}'.format(self.nw)))   
            
            name = 'sum_dimensions_{0}_{1}.txt'.format(self.begin,self.end-1)
            fmt = '%d ' + '%.4f '*len(self.sum_phase_dimension[0]) + '%.4f'
            dat = np.concatenate((np.array([range(self.begin,self.end)]).T, self.sum_phase_dimension, np.array([self.sum_rechi]).T),axis=1)
            np.savetxt(name,dat,fmt=fmt,
                        header='{0} \n {1}'.format(self.xyz,'id,  ls, ll, lw, lx,ly,lz, z0, z1_unpbc, chi'))   
        
        if self.is_analyze_proximity:
            name = 'sum_proximity_{0}_{1}.txt'.format(self.begin,self.end-1)
            dat = np.concatenate((np.array([range(self.begin,self.end)]).T, self.sum_prox),axis=1)
            fmt = '%d ' + '%d %d %d    '*len(self.ele_atomic_number) +'    %.4f'*2 # for each element, we have three phases
            np.savetxt(name,dat,fmt=fmt,
                        header='{0} \n {3} \n {1} \n {2}'.format(self.xyz,'solid liquid interface','id '+'           '.join(self.ele_chemical_symbol) +'          lw      chi', 'nw = {}'.format(self.nw)))   
            
    def analyze(self):
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


        for fn in range(self.begin,self.end,self.step):
            sol_liq_inter, ls, ll, lw, lx,ly,lz, z0, z1_unpbc, redchi = self.analyze_single(fn)
            self.sum_elements_counts.append(sol_liq_inter)
            self.sum_phase_dimension.append([ls, ll, lw, lx,ly,lz, z0, z1_unpbc])
            self.sum_rechi.append(redchi)
            if self.is_analyze_proximity:
                self.sum_prox.append(self.out)
            if self.is_analyze_mol_frac:
                self.analyze_mol_frac_single()
        
        self.plot_mol_frac()

    def analyze_mol_frac_single(self):
        inter = self.inter
        density_field_i = cal_mol_frac(inter)
        mol_frac         = []
        ele_list         = inter.unique_ele
        for i in range(len(ele_list)):
            ele = ele_list[i]
            selected_rho, selected_indice = cal_local_rho(inter, density_field_i[ele])
            selected_prox                 = np.array(self.proximity)[selected_indice]
            mean_prox_i, mean_rho_i       = _average_prox_vs_rho(selected_prox,selected_rho,int(self.inter.ngrid[0]))
            mol_frac.append(mean_prox_i); mol_frac.append(mean_rho_i)
        if self.mol_frac_frame_count == 0:
            self.mol_frac = np.array(mol_frac).T

        else:
            # merge two mol_frac, the old one and new one, 
            self.mol_frac = np.concatenate((self.mol_frac,np.array(mol_frac).T),axis=0)
        self.mol_frac_frame_count += 1
    
    def plot_mol_frac(self):
        self.load_default_setting()
        horizontal_scale = len(self.inter.unique_ele)*3
        fig,ax = plt.subplots(len(self.inter.unique_ele),1,figsize=(6,horizontal_scale),sharex=True,sharey=False)
        plt.xlabel('Proximity' + r'$(\mathrm{\AA})$')
        ele_list         = self.inter.unique_ele

        for i in range(len(self.inter.unique_ele)):
            ele = ele_list[i]
            # ith element takes 2*i 2*i+1 col
            mean_prox_i = self.mol_frac[:,2*i]
            mean_rho_i  = self.mol_frac[:,2*i+1]
            ax[i].plot(mean_prox_i, mean_rho_i,'r.', label=ele,alpha=0.3)
            ax[i].legend();ax[i].grid()
            ax[i].set_ylabel('Atomic fraction (%)')
        plt.savefig('atomic_fraction_vs_proximity.png',dpi=300,bbox_inches='tight')
        np.savetxt('atomic_fraction.txt',np.array(self.mol_frac), header='{0} \n {1}'.format(self.xyz,'           '.join(ele_list) ))  
        
    def analyze_single(self,idx,solid_center0=-1):
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

        sol_liq_inter = []
        # build interface, mda_xyz does not contain cell dim info, ase_xyz cannot read by pytim, we have to combine both
        self.inter, self.k,ase_a,mda_a = cal_inter(self.ase_xyz,self.mda_xyz,idx,self.mesh, self.alpha, self.level,self.mode)
        
        # element and mass list in the corresponding order
        self.inter.unique_mass = _unique_unsorted(self.inter.all_atoms.masses)
        self.inter.unique_ele  = _unique_unsorted(self.inter.all_atoms.elements)
        
        # build density map and project to target direction
        self.rho_zi,self.zi,rho_av     = temporal_mass_density_proj_coarse(self.inter,project_axis=self.project_axis,plot=False,level=self.level)        
        # find solid_center, move the 'center' to the center
        if solid_center0 <0: # if solid center is not specified, use the default value
            solid_center0 = self.zi[np.argmax(self.rho_zi)]
        for  solid_center in [solid_center0, solid_center0 -2, solid_center0+2, solid_center0-1, solid_center0+1, solid_center0-3, solid_center0+3]:
            self.zi_new,self.rho_zi_new    = pbc_x(self.zi,self.rho_zi,solid_center)
            # fit to double sided GDS function
            self.result,z0,z1,w1 = fit_gds_double_sided(self.zi_new,self.rho_zi_new,plot=False, verbose=self.verbose)
            # print("     result.redchi",result.redchi)
            if self.result.redchi<0.2:
                ## store the best fitting params to params
                params = self.result.params
                break
            
        # if params is None:  ## changing solid center is not the key, let us adjust the other initital values
        #     pass
            
        # if assert_chi:
        #     assert result.redchi<0.5        
        # choose inteface to be w 
        w = self.nw*w1  
        # ensure z0 is close to solid center
        # if (z0 - solid_center) <= (z1 - solid_center):
        #     pass
        # else:
        #     tmp = z0
        #     z0 = z1
        #     z1 = tmp
        dim      = ase_a.get_cell_lengths_and_angles()[self.project_axis]

        if (z0>dim) and (z1>dim):  ## corner case --  when solid center happens to be super close to the edge
            z0 = z0 - dim
            z1 = z1 - dim
        elif (z0<0) and (z1>0):
            z0 = z0 + dim
            z1 = z1 + dim    
        elif (z0 - solid_center) <= (z1 - solid_center):
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
        for ii in range(len(self.ele_atomic_number)):
            num = self.ele_atomic_number[ii]
            indices = self.ele_indices[ii]
            ez = ase_a.positions[indices][:,self.project_axis]
            ez[ez<0]   += ez[ez<0] + dim
            ez[ez>dim] -= ez[ez>dim] - dim
            sol_ez, liq_ez = count_z(z0,z1_unpbc, ez,w,dim)
            sol_liq_inter.append(len(sol_ez))
            sol_liq_inter.append(len(liq_ez))
            sol_liq_inter.append(len(indices) - len(sol_ez)-len(liq_ez))

            # sol_liq_inter.append([len(sol_ez), len(liq_ez), len(indices) - len(sol_ez)-len(liq_ez)])

        if z0 <  z1_unpbc:
            ls,ll = dim -(z1_unpbc - z0) - w, (z1_unpbc - z0) - w
        else:
            ls, ll = z0-z1_unpbc - w, dim - (z0-z1_unpbc) - w
            
        self.lw = self.nw*w
        
        if ls <0:
            ls = 0.00001
        if ll <0:
            ll = 0.00001
        lx,ly,lz = ase_a.get_cell_lengths_and_angles()[:3]
        self.solid_center = solid_center
        self.z0 = z0
        self.z1_unpbc = z1_unpbc
        self.w = w
        
        if self.is_analyze_proximity:
            self.proximity = cal_proximity(self.inter)[0]
            selected_rho, selected_indice =cal_local_rho(self.inter)
            selected_prox = np.array(self.proximity)[selected_indice]
            self.mean_prox,self.mean_rho = _average_prox_vs_rho(selected_prox,selected_rho,self.inter.ngrid[0])
            result, z0, w = fit_gds_single_sided(self.mean_prox,self.mean_rho,vary_z0=False,plot=False,verbose=False)
            
            self.proximity_lw = self.nw*w  # this w is different
            
            phase1 = self.inter.all_atoms.elements[tuple([self.proximity <= -self.proximity_lw/2])]
            phase2 = self.inter.all_atoms.elements[tuple([self.proximity >=self.proximity_lw/2])]
            interface = self.inter.all_atoms.elements[tuple([(self.proximity < self.proximity_lw/2) & (self.proximity > -self.proximity_lw/2)])]

            self.out = []
            for ele in self.ele_chemical_symbol:
                # out.append([len(phase2[phase2==ele]), len(phase1[phase1==ele]),  len(interface[interface==ele])])    
                self.out.append(len(phase2[phase2==ele]))
                self.out.append(len(phase1[phase1==ele]))
                self.out.append(len(interface[interface==ele]))
            self.out.append(self.proximity_lw)
            self.out.append(result.redchi)
                            
        return sol_liq_inter, ls, ll, self.lw, lx,ly,lz, self.z0, z1_unpbc, self.result.redchi

    def save_prox_xyz(self,outname='prox.xyz'):
        """
        save xyz file with k value being proximity or normalized 1,2,0

        Parameters
        ----------
        out : TYPE, optional
            DESCRIPTION. The default is 'prox.xyz'.

        Returns
        -------
        None.

        """
        out = open(outname,'w')
        
        
        cell = np.insert(self.inter.box, [1,1,1,2,2,2],[0,0,0,0,0,0])
        string_with_tabs = 'Lattice=\"' +" ".join(str(item) for item in cell) +'\"'
        ct = string_with_tabs + " Properties=species:S:1:pos:R:3:k:R:1\n" # k is the similarity kernel
        
        out.writelines(str(self.inter.all_atoms.n_atoms)+'\n')
        out.writelines(ct)
        
        lines = [' '.join(map(str, row)) + '\n' for row in self.inter.all_atoms.positions]
        
        phase1_array = np.array([int(item) for item in (self.proximity < -self.proximity_lw/2)])   # phase 1 k value is 1
        phase2_array = np.array([int(item) for item in (self.proximity >  self.proximity_lw/2)])*2 # phase 2 k value is 2
        phase_array = phase1_array + phase2_array  # interface is 0
        
        combined_lines = [self.inter.all_atoms.elements[i] + ' ' + line.strip() + ' ' + str(phase_array[i]) + '\n' for i, line in enumerate(lines)]
        out.writelines(combined_lines)
        out.close()
        
        # save raw proximity as k value
        out = open('raw_'+outname,'w')
        out.writelines(str(self.inter.all_atoms.n_atoms)+'\n')
        out.writelines(ct)
        combined_lines = [self.inter.all_atoms.elements[i] + ' ' + line.strip() + ' ' + str(self.proximity[i]) + '\n' for i, line in enumerate(lines)]
        out.writelines(combined_lines)
        out.close()
        
    def load_default_setting(self,fontsize=13,font_family='Arial'):
        """
        """
        import matplotlib as mpl
        mpl.rc('font',family = font_family)
        mpl.rc('font',size = fontsize)
        print("set font family as {0}".format(font_family))
        print("set font size as {0}".format(fontsize))
    
    def plot_prox(self,out='prox.png'):
        """
        

        Parameters
        ----------
        out : TYPE, optional
            DESCRIPTION. The default is 'prox.png'.

        Returns
        -------
        None.

        """
        self.load_default_setting()
        import matplotlib.pyplot as plt
        result, z0, w = fit_gds_single_sided(self.mean_prox,self.mean_rho,vary_z0=False,plot=False,verbose=True)

        # plt.plot(self.mean_prox,self.mean_rho)# use clip_on=False if need marker atop axes
        xrange = [min(self.mean_prox), max(self.mean_prox)]
        yrange = [min(self.mean_rho), max(self.mean_rho)]
        plt.figure()
        plt.plot(self.mean_prox,self.mean_rho,'o')
        plt.plot(self.mean_prox,result.best_fit,'-')
        plt.fill_between([z0-self.proximity_lw/2,z0+self.proximity_lw/2], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        plt.ylim(yrange)
        plt.xlim(xrange)
        plt.xlabel('Proximity' + r'$(\mathrm{\AA})$')
        plt.ylabel(r'$\rho$'+' ' + r'$(\mathrm{g}/\mathrm{cm}^{3})$')
        plt.minorticks_on()
        plt.savefig(out,dpi=300,bbox_inches='tight')
    
    def plot_gds(self,out='gds.png'):
        """
        Show the density profile, fitted results, and GDS location
        """
        import matplotlib.pyplot as plt
        self.load_default_setting()
        yrange = [min(self.rho_zi)-.5, max(self.rho_zi)+.5]
        xrange = [0,max(self.zi)]
        dim = xrange[1]
        z0s = [self.z0,self.z0]
        z1s = [self.z1_unpbc,self.z1_unpbc]       
        plt.figure()
        plt.plot(self.zi, self.rho_zi,'bo')
        plt.plot(z0s,yrange,'k--')
        plt.plot(z0s+dim,yrange,'k--')
        plt.plot(z0s-dim,yrange,'k--')

        plt.plot(z1s,yrange,'k--')
        plt.plot(z1s+dim,yrange,'k--')
        plt.plot(z1s-dim,yrange,'k--')

        rho_zi_new = self.result.best_fit
        # more flexibility
        # rho_zi_new = self.result.eval(self.result.params,z=self.zi_new)
        plt.plot(self.zi_new, rho_zi_new,'k-')
        plt.plot(self.zi_new+dim,rho_zi_new,'k-')
        plt.plot(self.zi_new-dim, rho_zi_new,'k-')        
        plt.plot([self.zi_new[-1]-dim,self.zi_new[0]],[rho_zi_new[-1], rho_zi_new[0]], 'k-')
    
        plt.fill_between([self.z0-self.lw/2,self.z0+self.lw/2], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        plt.fill_between([self.z1_unpbc-self.lw/2,self.z1_unpbc+self.lw/2], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        # consider PBC
        plt.fill_between([self.z0-self.lw/2 - dim,self.z0+self.lw/2 -dim], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        plt.fill_between([self.z0-self.lw/2 + dim,self.z0+self.lw/2 +dim], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
    
        plt.fill_between([self.z1_unpbc-self.lw/2-dim,self.z1_unpbc+self.lw/2-dim], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        plt.fill_between([self.z1_unpbc-self.lw/2+ dim,self.z1_unpbc+self.lw/2+ dim], [yrange[0], yrange[0]],[yrange[1], yrange[1]],color='r',alpha=0.3)
        plt.ylim(yrange)
        plt.xlim(xrange)
        plt.minorticks_on()
        plt.xlabel('z ' + r'$(\mathrm{\AA})$')
        plt.ylabel(r'$\rho$'+' ' + r'$(\mathrm{g}/\mathrm{cm}^{3})$')
        plt.savefig(out,dpi=300,bbox_inches='tight')
