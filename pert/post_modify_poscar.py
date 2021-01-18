#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 20:23:32 2021
use 
merge_out.py
extract_deepmd

1- dp test with model prefix
2- summarize dp test results compare energy, force, and virial

example:
    python ~/script/mldp/check_nbands_nelm.py -ip all
    python ~/script/mldp/merge_out.py -o OUTCAR -r y
    python ~/script/mldp/extract_deepmd.py -bs 1000 -o OUTCAR -d deepmd -t 
(deepmd-kit-gpu)[jd848@n7214 deepmd]$ dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/mpmt/mm12/pv_Jan8_cpu.pb -d mm12

dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/mpmt/mm12/pv_Jan8_cpu.pb -d mm12
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/extreme_filtered/re4/pv-cpu.pb -d extreme4  
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/extreme_filtered/re7/pv_gpu.pb -d extreme7   
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/pert/mp1/pv.pb -d mp1 
dp test -m /u/project/ESS/lstixrud/jd848/pv_hf_copy/dp-train/onlymelt/om2/pv.pb -d om2


@author: jiedeng
"""

import argparse
parser = argparse.ArgumentParser("Test model on poscars where atoms are manually palced and interatomic distances are controlled ")
parser.add_argument("--index","-i",nargs="+", help="index of atom of interest, if not provided, read from log file")
parser.add_argument("--log","-log", help="log file")
parser.add_argument("--test_folder","-tf",default='deepmd', help="where is dp test located?")
parser.add_argument("--model_prefix","-mp",nargs="+", help="model prefix, can take multiple ones")
parser.add_argument("--xrange","-xr",nargs="+",type=float, help="x range")
parser.add_argument("--yrange","-yr",nargs="+",type=float, help="y range")

args = parser.parse_args()

import numpy as np
import os
import matplotlib.pyplot as plt

if args.index:
    atom1, atom2 = args.index[0], args.index[1]
else:
    if args.log:
        tmp = np.loadtxt(args.log).astype(int)
    else:
        tmp = np.loadtxt('log').astype(int)
    atom1, atom2 = tmp[0], tmp[1]

def _make_dir(test_folder,prefix):
    e_dir = os.path.join(test_folder,prefix+'.e.out')
    f_dir = os.path.join(test_folder,prefix+'.f.out')
    v_dir = os.path.join(test_folder,prefix+'.v.out')

    e=np.loadtxt(e_dir)
    f=np.loadtxt(f_dir)
    v=np.loadtxt(v_dir)
    
    nframes = len(e)
    natoms  = len(f)//nframes
    
    coord_dir = os.path.join(os.path.join(test_folder,'set.000'),'coord.npy')
    box_dir   = os.path.join(os.path.join(test_folder,'set.000'),'box.npy')
    coord     = np.load(coord_dir)
    box       = np.load(box_dir)[0]
    box       = np.array([box[0], box[4], box[8]])
    return e, f, v, coord, box, nframes, natoms

def _interatomic_distance(coord, box, atom1, atom2):
    coord_atom1 = coord[:,(atom1*3):(atom1*3+3)]
    coord_atom2 = coord[:,(atom2*3):(atom2*3+3)]
    
    #cal distance 
    diffs = coord_atom2 - coord_atom1
    cond = np.where(diffs > box / 2.)
    diffs[cond] -= box[cond[1]]
    cond = np.where(diffs < -box / 2.)
    diffs[cond] += box[cond[1]]
    sqdist = np.square(diffs).sum(axis=1)
    dist = np.sqrt(sqdist)   
    return dist

def _extract_force(f,atom1, atom2, nframes, natoms):
    f_idx = np.array(range(nframes))*natoms + atom2
    f2 = f[f_idx,:] #force felt by atom2
    f2_vasp = f2[:,:3]; mag_f2_vasp = np.sqrt(np.sum(f2_vasp**2,axis=1))
    f2_nn   = f2[:,3:]; mag_f2_nn   = np.sqrt(np.sum(f2_nn**2,axis=1))
    return mag_f2_vasp, mag_f2_nn

def _extract_e(e):
    return e[:,0], e[:,1]

def get_dist_f_e(test_folder, prefix, atom1, atom2):
    e, f, v, coord, box, nframes, natoms = _make_dir(test_folder,prefix)
    dist = _interatomic_distance(coord, box, atom1, atom2)
    mag_f2_vasp, mag_f2_nn = _extract_force(f,atom1, atom2, nframes, natoms)
    e_vasp, e_nn = _extract_e(e)
    return dist, mag_f2_vasp, mag_f2_nn, e_vasp, e_nn

def autoscale_y(ax,margin=0.1):
    """
    This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims
    expample:
        import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(-100,100,1000)
    y = x**2 + np.cos(x)*100

    fig,axs = plt.subplots(1,2,figsize=(8,5))

    for ax in axs:
        ax.plot(x,y)
        ax.plot(x,y*2)
        ax.plot(x,y*10)
        ax.set_xlim(-10,10)

    autoscale_y(axs[1])

    axs[0].set_title('Rescaled x-axis')
    axs[1].set_title('Rescaled x-axis\nand used "autoscale_y"')

    plt.show()
    """


    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)
    

count = 0
fig,ax = plt.subplots(2,1,figsize=(6,8),sharex=True)
for prefix in args.model_prefix:

    dist, mag_f2_vasp, mag_f2_nn, e_vasp, e_nn = \
    get_dist_f_e(args.test_folder, prefix, atom1, atom2)

    if count == 0:
        ax[0].plot(dist,mag_f2_vasp,'o',label='vasp')
        ax[1].plot(dist,e_vasp,'o',label='vasp')       
    count += 1
    
    ax[0].plot(dist,mag_f2_nn,'+',label=prefix)
    ax[1].plot(dist,e_nn,'+',label=prefix)

ax[0].set_ylabel('force (eV/A)')
ax[1].set_ylabel('energy (eV)')
ax[1].set_xlabel('interatomic distance (A)')
ax[0].grid()
ax[1].grid()
ax[0].legend()
ax[1].legend()
if args.xrange:
    plt.xlim(args.xrange)
    autoscale_y(ax[0])
    autoscale_y(ax[1])
plt.minorticks_on()
plt.show()    


