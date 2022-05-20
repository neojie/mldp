#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:26:56 2022

@author: jiedeng
"""
import numpy as np
log = open('log.sub','r')
count = 0
for line in log:
    if count ==0:
        line0 = line
        ch = np.loadtxt(line.strip().replace('sh','txt'))
    else:
        ch = np.concatenate((ch,np.loadtxt(line.strip().replace('sh','txt'))))
    count += 1
    
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

import os
name='stat_{0}_{1}.txt'.format(int(ch[:,0][0]),int(ch[:,0][-1]))
save_ch(ch,0,2.5,os.path.abspath('.'),name)