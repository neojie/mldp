# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
""" Module: pytim.gaussian_kde_pbc
    ==============================
"""
from __future__ import print_function
import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial import cKDTree


class gaussian_kde_pbc(gaussian_kde):
    # note that here "points" are those on the grid

    def evaluate_pbc_fast(self, points,atomic_masses):
        grid = points
        
        pos = self.pos
        box = self.box
        d = self.sigma * 2.5
#        print("d is",d)
#        d = self.sigma
        results = np.zeros(grid.shape[1], dtype=float)
        results2 = np.zeros(grid.shape[1], dtype=float)
        gridT = grid.T[:]
#        print("gridT shape is",np.shape(gridT))
#        print("gridT is",gridT)
        tree = cKDTree(gridT, boxsize=box)
#        tree = cKDTree(gridT)

        # the indices of grid elements within distane d from each of the pos
        scale = 2. * self.sigma**2
        indlist = tree.query_ball_point(pos, d)
#        print("shape of indlist",indlist.shape)
        for n, ind in enumerate(indlist):
            dr = gridT[ind, :] - pos[n]
            #best explanation for np.where https://stackoverflow.com/questions/34667282/numpy-where-detailed-step-by-step-explanation-examples
            cond = np.where(dr > box / 2.)
            dr[cond] -= box[cond[1]]
            cond = np.where(dr < -box / 2.)
            dr[cond] += box[cond[1]]
            #make sure dr not exceed half the box
            dens      = np.exp(-np.sum(dr * dr, axis=1) / scale)
            mass_dens = np.exp(-np.sum(dr * dr, axis=1) / scale)*atomic_masses[n]
            results[ind] += dens
            results2[ind] += mass_dens
        NA      = 6.022e23
        volume_tot  = self.box[0]*self.box[1]*self.box[2]*1e-24 # cm^3
        volume0     = volume_tot/len(results) # cm^3/mol
        norm_factor = sum(atomic_masses)/NA/sum(results2*volume0)
        results2 = results2*norm_factor
        
        return results,results2

    
    def evaluate_mass_pbc_fast(self, points):
        """
        added by Jie
        """
        grid = points
        pos = self.pos
        box = self.box
        d = self.sigma * 2.5
#        print("d is",d)
#        d = self.sigma

        results = np.zeros(grid.shape[1], dtype=float)
        gridT = grid.T[:]
#        print("gridT shape is",np.shape(gridT))
#        print("gridT is",gridT)
        tree = cKDTree(gridT, boxsize=box)
#        tree = cKDTree(gridT)

        # the indices of grid elements within distane d from each of the pos
        scale = 2. * self.sigma**2
        indlist = tree.query_ball_point(pos, d)
#        print("shape of indlist",indlist.shape)
        for n, ind in enumerate(indlist):
            dr = gridT[ind, :] - pos[n]
            #best explanation for np.where https://stackoverflow.com/questions/34667282/numpy-where-detailed-step-by-step-explanation-examples
            cond = np.where(dr > box / 2.)
            dr[cond] -= box[cond[1]]
            cond = np.where(dr < -box / 2.)
            dr[cond] += box[cond[1]]
            #make sure dr not exceed half the box
            dens = np.exp(-np.sum(dr * dr, axis=1) / scale)*self.universe.atoms.masses[n]
            results[ind] += dens

        return results
