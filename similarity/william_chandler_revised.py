#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:14:43 2024

@author: jiedeng
"""

from pytim import WillardChandler
from pytim import utilities
from pytim.sanity_check import SanityCheck
from pytim import messages
from pytim.patches import patchTrajectory, patchOpenMM, patchMDTRAJ
from pytim.vtk import Writevtk
from skimage import measure
import numpy as np
from gaussian_kde_pbc import gaussian_kde_pbc

try:
    marching_cubes = measure.marching_cubes
except AttributeError:
    marching_cubes = measure.marching_cubes_lewiner
    
def density_map(pos, grid, sigma, box):
    values = np.vstack([pos[::, 0], pos[::, 1], pos[::, 2]])
    kernel = gaussian_kde_pbc(values, bw_method=sigma / values.std(ddof=1))
    kernel.box = box
    kernel.sigma = sigma
    return kernel, values.std(ddof=1)
   
class WillardChandlerRevised(WillardChandler):
    def __init__(self,
                 universe,
                 group=None,
                 alpha=2.0,
                 radii_dict=None,
                 mesh=2.0,
                 symmetry='spherical',
                 cluster_cut=None,
                 cluster_threshold_density=None,
                 extra_cluster_groups=None,
                 centered=False,
                 warnings=False,
                 autoassign=True,
                 level=None,
                 **kargs):
        self.level = level
        self.autoassign, self.do_center = autoassign, centered
        sanity = SanityCheck(self, warnings=warnings)
        sanity.assign_universe(universe, group)
        sanity.assign_alpha(alpha)

        if mesh <= 0:
            raise ValueError(messages.MESH_NEGATIVE)
        self.mesh, self.spacing, self.ngrid, self.PDB = mesh, None, None, {}

        sanity.assign_radii(radii_dict=radii_dict)

        sanity.assign_cluster_params(cluster_cut,
                                      cluster_threshold_density, extra_cluster_groups)

        self._assign_symmetry(symmetry)
        patchTrajectory(self.universe.trajectory, self)
        self._assign_layers()
        self._atoms = self._layers[:]  # this is an empty AtomGroup
        self.writevtk = Writevtk(self)

        
    def _assign_layers(self):
        """ There are no layers in the Willard-Chandler method.

            This function identifies the dividing surface and stores the
            triangulated isosurface, the density and the particles.

        """
        self.reset_labels()
        # we assign an empty group for consistency
        self._layers, self.normal = self.universe.atoms[:0], None

        self.prepare_box()

        self._define_cluster_group()

        self.centered_positions = None
        if self.do_center is True:
            self.center()

        pos = self.cluster_group.positions
        box = self.universe.dimensions[:3]

        ngrid, spacing = utilities.compute_compatible_mesh_params(
            self.mesh, box)
        self.spacing, self.ngrid = spacing, ngrid
        grid = utilities.generate_grid_in_box(box, ngrid, order='xyz')
        self.kernel, _ = utilities.density_map(pos, grid, self.alpha, box)

        self.kernel.pos = pos.copy()
        self.grid = grid
        self.box = box
        #self.density_field = kernel.evaluate_pbc_fast(grid) #JD
        self.density_field,self.mass_density_field = self.kernel.evaluate_pbc_fast(self.grid,self.universe.atoms.masses) #JD
        # print("here,", self.mass_density_field)
        ## added by Jie, normalize mass density
        # NA      = 6.022e23
        # volume_tot  = self.box[0]*self.box[1]*self.box[2]*1e-24
        # volume0     = volume_tot/(self.ngrid[-1]*self.ngrid[0]*self.ngrid[1]) # cm^3/mol
        # norm_factor = sum(self.universe.atoms.masses)/NA/sum(self.mass_density_field*volume0)
        # self.mass_density_field = self.mass_density_field*norm_factor
        # print("here1,", self.mass_density_field)
        
        # Thomas Lewiner, Helio Lopes, Antonio Wilson Vieira and Geovan
        # Tavares. Efficient implementation of Marching Cubes’ cases with
        # topological guarantees. Journal of Graphics Tools 8(2) pp. 1-15
        # (december 2003). DOI: 10.1080/10867651.2003.10487582
        # volume = self.density_field.reshape(
        #     tuple(np.array(ngrid[::-1]).astype(int)))
        # verts, faces, normals, values = marching_cubes(
        #     volume, None, spacing=tuple(spacing))
        self.mass_density_field_reshape = self.mass_density_field.reshape(
            tuple(np.array(ngrid).astype(int)))
        verts, faces, normals, values = marching_cubes(
            self.mass_density_field_reshape, self.level, spacing=tuple(spacing))
        # note that len(normals) == len(verts): they are normals
        # at the vertices, and not normals of the faces
        # verts and normals have x and z flipped because skimage uses zyx ordering
        # self.triangulated_surface = [np.fliplr(verts), faces, np.fliplr(normals)]
        self.triangulated_surface = [verts, faces, normals]

        self.surface_area = measure.mesh_surface_area(verts, faces)
        verts += spacing[::-1] / 2.
