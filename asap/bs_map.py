#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 20:22:07 2020

@author: jiedeng
"""
#from asaplib.reducedim import KernelPCA





""" for maps and fits """

"""
/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/ASAP/build/lib/asaplib/cli/func_asap.py
"""
def read_xyz_n_dm(fxyz, design_matrix, use_atomic_descriptors, only_use_species, peratom):
    dm = []
    dm_atomic = []
    # try to read the xyz file
    if fxyz is not None and fxyz != 'none':
        from asaplib.data import ASAPXYZ
#        from asapxyzs import ASAPXYZs
        asapxyz = ASAPXYZ(fxyz)
        if use_atomic_descriptors:
            dm = asapxyz.get_atomic_descriptors(design_matrix, only_use_species)
        else:
            dm, dm_atomic = asapxyz.get_descriptors(design_matrix, peratom)
    else:
        asapxyz = None
        print("Did not provide the xyz file. We can only output descriptor matrix.")
    # we can also load the descriptor matrix from a standalone file
    import os
    if os.path.isfile(design_matrix[0]):
        try:
            import numpy as np
            dm = np.genfromtxt(design_matrix[0], dtype=float)
            print("loaded the descriptor matrix from file: ", design_matrix[0])
        except:
            raise ValueError('Cannot load the descriptor matrix from file')
    return asapxyz, dm, dm_atomic

"""for maps"""


def figure_style_setups(prefix,
                        colorlabel, colorscale, colormap,
                        style, aspect_ratio, adjusttext):
    fig_options = {'outfile': prefix,
                   'show': False,
                   'title': None,
                   'size': [8 * aspect_ratio, 8],
                   'cmap': colormap,
                   'components': {
                       'first_p': {'type': 'scatter', 'clabel': colorlabel,
                                   'vmin': colorscale[0], 'vmax': colorscale[1]},
                       'second_p': {"type": 'annotate', 'adtext': adjusttext}}
                   }
    if style == 'journal':
        fig_options.update({'xlabel': None, 'ylabel': None,
                            'xaxis': False, 'yaxis': False,
                            'remove_tick': True,
                            'rasterized': True,
                            'fontsize': 12,
                            'size': [4 * aspect_ratio, 4]
                            })
    return fig_options


def map_process(obj, reduce_dict, axes, map_name):
    """
    process the dimensionality reduction command
    """
    # project
    if 'type' in reduce_dict.keys() and reduce_dict['type'] == 'RAW':
        proj = obj['design_matrix']
        if obj['map_options']['peratom']:
            proj_atomic = obj['design_matrix_atomic']
        else:
            proj_atomic = None
    else:
        from asaplib.reducedim import Dimension_Reducers
        dreducer = Dimension_Reducers(reduce_dict)
        proj = dreducer.fit_transform(obj['design_matrix'])
        if obj['map_options']['peratom']:
            print("Project atomic design matrix with No. of samples:", len(obj['design_matrix_atomic']))
            proj_atomic = dreducer.transform(obj['design_matrix_atomic'])
        else:
            proj_atomic = None
    # plot
    fig_spec = obj['fig_options']
    plotcolor = obj['map_options']['color']
    plotcolor_atomic = obj['map_options']['color_atomic']
    annotate = obj['map_options']['annotate']
    if 'cluster_labels' in obj.keys():
        labels = obj['cluster_labels']
    else:
        labels = []
    map_plot(fig_spec, proj, proj_atomic, plotcolor, plotcolor_atomic, labels, annotate, axes)
    # output 
    outfilename = obj['fig_options']['outfile']
    outmode = obj['map_options']['outmode']
    species_name = obj['map_options']['only_use_species']
    if obj['map_options']['project_atomic']:
        map_save(outfilename, outmode, obj['asapxyz'], None, proj, map_name, species_name)
    else:
        map_save(outfilename, outmode, obj['asapxyz'], proj, proj_atomic, map_name, species_name)


def map_plot(fig_spec, proj, proj_atomic, plotcolor, plotcolor_atomic, labels, annotate, axes):
    """
    Make plots
    """
    from matplotlib import pyplot as plt
    from asaplib.plot import Plotters
    asap_plot = Plotters(fig_spec)
    asap_plot.plot(proj[:, axes], plotcolor[:], labels[:], annotate[:])
    if proj_atomic is not None:
        asap_plot.plot(proj_atomic[:, axes], plotcolor_atomic[:], [], [])
    plt.show()


def map_save(foutput, outmode, asapxyz, proj, proj_atomic, map_name, species_name):
    """
    Save the low-D projections
    """
    if outmode == 'matrix':
        import numpy as np
        if proj is not None:
            np.savetxt(foutput + ".coord", proj, fmt='%4.8f', header='low D coordinates of samples')
        if proj_atomic is not None:
            np.savetxt(foutput + "-atomic.coord", proj_atomic, fmt='%4.8f', header=map_name)
    elif outmode in ('xyz', 'chemiscope'):
        if proj is not None:
            asapxyz.set_descriptors(proj, map_name)
        if proj_atomic is not None:
            asapxyz.set_atomic_descriptors(proj_atomic, map_name, species_name)
        if outmode == 'xyz':
            asapxyz.write(foutput)
        else:
            # If we write atomic projection assume we want to show them
            cutoff = 3.5 if proj_atomic else None
            asapxyz.write_chemiscope(foutput, cutoff=cutoff)
    else:
        pass
    

"""
/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/ASAP/build/lib/asaplib/cli/cmd_asap.py
"""




"""
/Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/ASAP/asaplib/plot/plot_colors.py
"""
import os
import numpy as np
def set_color_function(fcolor='none', asapxyz=None, colorscol=0, n_samples=0, 
              peratom=False, project_atomic=False, use_atomic_species=None, color_from_zero=False, extensive=False):
    """ obtain the essential informations to define the colors of data points
    Parameters
    ----------
    fcolor: str
             the name of the file or the tag in the xyz to define the colors
    asapxyz: ASAPXYZ object, (optional)
    colorscol: int, (optional). 
              if the color file has more than one column, which column to use
    n_samples: int, (optional). 
              The number of data points
    peratom: bool
              return atomic color
    project_atomic: bool
              the samples are atomic descriptors
    use_atomic_species: int
              the atomic number of the selected species
    color_from_zero: bool
              set the min color to zero
    extensive: bool
              normalize the quatity by number of atoms
    """

    plotcolor = []
    plotcolor_atomic = []
    colorscale = [None, None]

    # if there is a file named "fcolor", we load it for the color scheme
    if os.path.isfile(fcolor):
        # load the column=colorscol for color functions
        try:
            loadcolor = np.genfromtxt(fcolor, dtype=float)
        except:
            raise IOError('Error in loading fcolor files for the color scheme')

        # print(np.shape(loadcolor))
        if colorscol > 0 or len(np.shape(loadcolor)) > 1:
            plotcolor = loadcolor[:, colorscol]
        else:
            plotcolor = loadcolor
        print('load file: ' + fcolor + ' for color schemes')

        if peratom or project_atomic:
            if asapxyz is None:
                raise IOError('Need the xyz so that we know the number of atoms in each frame')
            elif asapxyz.get_num_frames() == len(plotcolor):
                for index, natomnow in enumerate(asapxyz.get_natom_list_by_species(use_atomic_species)):
                    plotcolor_atomic = np.append(plotcolor_atomic, plotcolor[index] * np.ones(natomnow))
            elif asapxyz.get_total_natoms() == len(plotcolor):
                plotcolor_atomic = plotcolor
            else:
                raise ValueError('Length of the xyz trajectory is not the same as number of colors in the fcolor file')

    elif n_samples > 0 and (fcolor == None or fcolor == 'none' or fcolor == 'Index' or fcolor == 'index') and peratom == False:
        # we use the index as the color scheme
        plotcolor = np.arange(n_samples)
        fcolor = 'sample index'

    elif asapxyz is None:
        raise IOError('Cannot find the xyz or fcolor files for the color scheme')

    else:
        if fcolor == None or fcolor == 'none' or fcolor == 'Index' or fcolor == 'index':
            # we use the index as the color scheme
            plotcolor = np.arange(asapxyz.get_num_frames())
            fcolor = 'sample index'
            if peratom or project_atomic:
                for index, natomnow in enumerate(asapxyz.get_natom_list_by_species(use_atomic_species)):
                    plotcolor_atomic = np.append(plotcolor_atomic, plotcolor[index] * np.ones(natomnow))
        else:
            try:
                plotcolor = asapxyz.get_property(fcolor, extensive)
            except:
                raise ValueError('Cannot find the specified property from the xyz file for the color scheme')
            if peratom or project_atomic:
                try:
                    plotcolor_atomic = asapxyz.get_atomic_property(fcolor, extensive, [], use_atomic_species)
                    #print(np.shape(plotcolor_atomic))
                except:
                    raise ValueError('Cannot find the specified atomic property from the xyz file for the color scheme')

    if color_from_zero:
        # set the min to zero
        plotcolor -= np.ones(len(plotcolor))*np.nanmin(plotcolor)
        plotcolor_atomic -= np.ones(len(plotcolor_atomic))*np.nanmin(plotcolor)

    colorlabel = str(fcolor)
    if peratom and not project_atomic:
        # print(np.shape(plotcolor_atomic))
        colorscale = [np.nanmin(plotcolor_atomic), np.nanmax(plotcolor_atomic)]
        return plotcolor, np.asarray(plotcolor_atomic), colorlabel, colorscale
    elif project_atomic:
        colorscale = [None, None]
        return np.asarray(plotcolor_atomic), [], colorlabel, colorscale
    else:
        colorscale = [None, None]
        return plotcolor, [], colorlabel, colorscale


class COLOR_PALETTE:
    def __init__(self, style=1):
        if style == 1:
            self.pal = ["#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#0000A6", "#63FFAC", "#B79762",
                        "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                        "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                        "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF",
                        "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
                        "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99",
                        "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1",
                        "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459",
                        "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375",
                        "#A3C8C9", "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F",
                        "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9",
                        "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79",
                        "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534",
                        "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#7A4900"]
        elif style == 2:
            self.pal = ["#30a2da", "#fc4f30", "#e5ae38", "#6d904f", "#8b8b8b", "#006FA6", "#A30059", "#af8dc3",
                        "#922329", "#1E6E00"]

        self.n_color = len(self.pal)

    def __getitem__(self, arg):  # color cycler
        assert arg > -1, "???"
        return self.pal[arg % self.n_color]
    

def map(obj, fxyz, design_matrix, prefix, output, extra_properties,
        use_atomic_descriptors, only_use_species, peratom, keepraw,
        color, color_column, color_label, colormap, color_from_zero, normalized_by_size,
        annotate, adjusttext, style, aspect_ratio):
    """
    Making 2D maps using dimensionality reduction.
    This command function evaluated before the specific ones,
    we setup the general stuff here, such as read the files.
    """

    if not fxyz and not design_matrix[0]:
        return
    if prefix is None: prefix = "ASAP-lowD-map"
    obj['asapxyz'], obj['design_matrix'], obj['design_matrix_atomic'] = read_xyz_n_dm(fxyz, design_matrix,
                                                                                                  use_atomic_descriptors,
                                                                                                  only_use_species,
                                                                                                  peratom)

    # Read additional properties
    if extra_properties:
        obj['asapxyz'].load_properties(extra_properties)

    if obj['asapxyz'] is None: output = 'matrix'
    print(len(obj['design_matrix_atomic']))
    # remove the raw descriptors
    if not keepraw and obj['asapxyz'] is not None:
        print("Remove raw desciptors..")
        obj['asapxyz'].remove_descriptors(design_matrix)
        obj['asapxyz'].remove_atomic_descriptors(design_matrix)

    # color scheme
#    from asaplib.plot import set_color_function
    plotcolor, plotcolor_peratom, colorlabel, colorscale = set_color_function(color, obj['asapxyz'], color_column,
                                                                              0, peratom, use_atomic_descriptors,
                                                                              only_use_species, color_from_zero,
                                                                              normalized_by_size)
    if color_label is not None: colorlabel = color_label

    obj['map_options'] = {'color': plotcolor,
                              'color_atomic': plotcolor_peratom,
                              'project_atomic': use_atomic_descriptors,
                              'only_use_species': only_use_species,
                              'peratom': peratom,
                              'annotate': [],
                              'outmode': output,
                              'keepraw': keepraw
                              }
    if annotate != 'none':
        try:
            obj['map_options']['annotate'] = obj['asapxyz'].get_property(annotate)
        except:
            import numpy as np
            obj['map_options']['annotate'] = np.loadtxt(annotate, dtype="str")[:]

    obj['fig_options'] = figure_style_setups(prefix, colorlabel, colorscale, colormap, style, aspect_ratio,
                                                 adjusttext)
    return obj

def skpca(obj, scale, dimension, axes,
          kernel, kernel_parameter, sparse_mode, n_sparse):
    """Sparse Kernel Principal Component Analysis"""
    map_name = "skpca-d-" + str(dimension)
    reduce_dict = {}
    if scale:
        print("Perform standard scaling of the design matrix. To turn it off use `--no-scale`")
        reduce_dict = {"preprocessing": {"type": 'SCALE', 'parameter': None}}
    reduce_dict['skpca'] = {"type": 'SPARSE_KPCA',
                            'parameter': {"n_components": dimension,
                                          "sparse_mode": sparse_mode, "n_sparse": n_sparse,
                                          "kernel": {"first_kernel": {"type": kernel, "d": kernel_parameter}}}}
    map_process(obj, reduce_dict, axes, map_name)
    
## defaults based on /Users/jiedeng/opt/anaconda3/lib/python3.7/site-packages/ASAP/build/lib/asaplib/cli/cmd_cli_options.py
#obj = {}
##fxyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/recal/asap/ASAP-desc.xyz'
#fxyz = ['/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid0/r3/cont1/recal/asap/ASAP-desc.xyz', 
#        '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap/ASAP-desc.xyz']
##fxyz = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/asap/*'
##fxyz = ['/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/asap/k3r3recal.xyz', '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/asap/k1r1recal.xyz']
#design_matrix = '[*]'
#prefix = 'ASPA'
#output = 'none'
#extra_properties = False
#use_atomic_descriptors = False
#only_use_species = None
#peratom = False
#keepraw = False
#color = 'none'
#color_column = 0
#color_label = ''
#colormap = 'gnuplot'
#color_from_zero = False
#normalized_by_size = True
#annotate = 'none'
#adjusttext = False
#style = 'journal'
#aspect_ratio =2
#        
#obj1=map(obj, fxyz, design_matrix, prefix, output, extra_properties,
#        use_atomic_descriptors, only_use_species, peratom, keepraw,
#        color, color_column, color_label, colormap, color_from_zero, normalized_by_size,
#        annotate, adjusttext, style, aspect_ratio)
#
#
#scale =True#False benchmark with the command line 
#dimension = 5
#axes = [0, 1]
#kernel = 'linear'
#kernel_parameter = None
#sparse_mode = 'fps'
#n_sparse = 100
#
#skpca(obj1,scale, dimension, axes,
#          kernel, kernel_parameter, sparse_mode, n_sparse)