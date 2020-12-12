#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 09:41:11 2020


if target is 10 meV/atom, we so should cover every 10 meV/atom*Natom energy gap we have here is de
de/160/
@author: jiedeng
"""

import dpdata
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--outcar","-o",default='OUTCAR',help="OUTCAR")
parser.add_argument("--threshold","-t",type=float,default=2e-3,help="e/atom threshold")
parser.add_argument('--plot',"-p", default=True, action='store_false',help="plot the results?")

args   = parser.parse_args()

#outcar    = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/4k/100g/r3/OUTCAR'
#threshold = 2e-3 # 
#plot      = True
outcar = args.outcar
threshold = args.threshold
###############
ls  = dpdata.LabeledSystem(outcar,fmt='outcar')
e   = ls['energies']
de  = threshold*ls.get_natoms()
bins = np.arange(min(e),max(e+de),de)

idxs = []
for ii in range(len(bins)-1):
    lower = bins[ii]
    upper = bins[ii+1]    
    idx1  = np.where(e>=lower)[0]
    idx2  = np.where(e<upper)[0]
    idx   = np.intersect1d(idx1,idx2)
    try:
        idxs.append(idx[len(idx)//2])
#        idxs.append(idx[-1])
    except:
        pass
idxs=np.array(idxs).astype(int)
e_idxs=e[idxs]
assert len(np.unique(idxs)) == len(idxs) # idx should not contains duplicate values
print('{0}/{1} frames are selected with e/atom = {2}'.format(len(e_idxs),len(e),threshold))
idxs.sort()
idxfile = 'idx_unique_e'
np.savetxt(idxfile,idxs,fmt='%d',header='threshold={0}'.format(threshold))
print('idx file saved as {0}'.format(idxfile))


if args.plot:
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(1,2,figsize=(10,4))
    ax[0].hist(e,len(bins))
    ax[0].hist(e_idxs,len(bins))
    ax[0].set_xlabel('ETOT(eV)')
    ax[0].set_ylabel('count')
    ax[1].plot(e)
    ax[1].plot(idxs,e[idxs],'.')
    ax[1].set_ylabel('ETOT(eV)')
    ax[1].set_xlabel('step')
    fig.savefig('unique_e.png')   
    plt.show()

    


