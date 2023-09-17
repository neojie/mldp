#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:25:43 2020
## msd direct is wrong with unwrapping
@author: jiedeng
"""

import time
start_time = time.time()
formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
print('Program starts at: ', formatted_time)
print("Loading modules ...")
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import argparse
import glob
import multiprocessing
import os
module_time = time.time()
print("End after %.2f s" % (module_time - start_time))

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str,help="input file")
parser.add_argument("--eles","-e",type = str,default = None, help="elements to analyze")
parser.add_argument("--format","-ft",type = str,default = 'LAMMPSDUMP', help="file format, e.g., LAMMPSDUMP, PDB")
parser.add_argument("--timestep","-ts",type = int,default = 1, help="timestep")

args   = parser.parse_args()


#def unwrap(idx):
#    #Bulow et al., 2020 unwrap method eqn3
#    newcoords = []
#    for i in range(len(frames)):
#        if i == 0:
#            newcoords.append(frames[i].positions[idx])
#        else:
#            curr = frames[i].positions[idx]
#            
#            newcoords.append(curr - np.floor((curr - newcoords[i-1])/box +.5)*box)
#
#    return np.array(newcoords)


def unwrap(idx):
    #Bulow et al., 2020 unwrap method, eqn1
    newcoords = []
    for i in range(len(frames)):
        if i == 0:
            newcoords.append(frames[i].positions[idx])
        else:
            wr_i = frames[i].positions[idx]
            wr_i1 = frames[i-1].positions[idx]
            un_i1 = newcoords[i-1]
            tmp =   np.floor((wr_i - wr_i1)/box +.5)*box
            newcoords.append(un_i1 + wr_i - wr_i1 - tmp)
            
    return np.array(newcoords)


def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_fft(r):
  N=len(r)
  D=np.square(r).sum(axis=1) 
  D=np.append(D,0) 
  S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
  Q=2*D.sum()
  S1=np.zeros(N)
  for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
  return S1-2*S2

def msd_straight_forward(r):
    shifts = np.arange(len(r))
    msds = np.zeros(shifts.size)    

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()

    return msds

def task(i):
    r = unwrap(i)
    msd_tmp = msd_fft(r)
#   msd_tmp =  msd_straight_forward(r)
    print(i, end=' ')
    return msd_tmp

if args.file:
    file = args.file
else:
    print("No file supplied, search for *dump file")
    files  = glob.glob("*dump")
    file = files[0]
    print("Find {0}; analyze {1}".format(files,file))

#file = '/Users/jiedeng/Documents/ml/deepmd-kit/my_example/6k/tmp/r1/20/mgsio3.dump'
print("Loading dump file ...")
u_md = mda.Universe(file,format=args.format)
loading_time = time.time()
print("End after %.2f s" % (loading_time - module_time))

print("Transferring to memory ...")
u_md.transfer_to_memory()
transfer_time = time.time()
print("End after %.2f s" % (transfer_time - loading_time))

frames = u_md.trajectory
box    = frames[0].dimensions[0:3]

eles = []
for atom in u_md.atoms.types:
    if not (atom in eles):
        eles.append(atom)

# find unique elements
if not (args.eles):
    ele_sel = eles
else:
    ele_sel = args.eles.split('-')

print("Calculating MSD ...")
if cores := os.getenv('SLURM_CPUS_PER_TASK'):
    cores = int(cores)
else:
    cores = multiprocessing.cpu_count()
    if cores > 16:
        cores = 16
print('Number of CPU cores', cores)
msds = []
for ele in ele_sel:
    idx = u_md.select_atoms('type {0}'.format(ele)).indices
    print()
    print('Elements:', ele)
    print('Number of atoms:',len(idx))
    tmp = np.zeros(len(u_md.trajectory))
    with multiprocessing.Pool(cores) as pool:
        results = pool.map(task, idx)
    for result in results:
        tmp += result
    msds.append(tmp/len(idx))
print()
print(msds)
calc_time = time.time()
print("End after %.2f s" % (calc_time - transfer_time))

print("Saving to file ...")
np.savetxt("msd_fft.txt", np.array(msds).T, header='    '.join(ele_sel), fmt = '%2.6f')
sf_time = time.time()
print("End after %.2f s" % (sf_time - calc_time))

print("Plotting ...")
plt.figure()
for i in range(len(ele_sel)):
    plt.loglog((np.array(list(range(len(msds[i]))))*args.timestep)[1:],msds[i][1:],label=eles[i])
plt.legend()
#plt.ylim([1e-1,1e2])
plt.savefig("msd_fft.png",bbox_inches='tight')
plot_time = time.time()
print("End after %.2f s" % (plot_time - sf_time))
