#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 19:16:57 2020

@author: jiedeng
"""
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Cv for vasp")
parser.add_argument("--input_dir","-i", type=str,default='.', help="out directory, default: cwd")
parser.add_argument("--beg","-b", type=int,default=0, help="begin from index")
parser.add_argument("--num","-n", type=int, help="number of atoms in the system")

args = parser.parse_args()

input_dir = args.input_dir
outcar   = os.path.join(input_dir,'OUTCAR')
fo      = open(outcar)
lines   = fo.readlines()


get_eentropy = False
get_temperature = False
get_e = False
get_v = False
vol = 0 
dat = []

for i in range(1000,len(lines)):  # can start from a large number since typically OUTCAR first 5000 lines does not contatin that info
    line = lines[i]
    if not get_v:
        if 'volume of cell' in line:
            vol = float(line.split()[-1])
            get_v = True
    if 'EENTRO' in line:
        get_eentropy = True
        eentropy     = float(line.split()[-1])
    if 'mean temperature' in line:
        get_temperature = True
        temperature =  float(line.split()[-1])
    if 'energy without entropy' in line:
        get_e = True
        e = float(line.split()[4])
    if get_eentropy and get_temperature and get_e:
        dat.append([temperature,eentropy,e])
        get_eentropy, get_temperature,get_e = False, False,False

fo.close()

outcar   = os.path.join(input_dir,'OUTCAR')
fo      = open(outcar)
lines   = fo.readlines()

dat = np.array(dat)[args.beg:,:]

#dat = np.array(dat)[3000:,:]


def fluct(data):
    tmp=data - data.mean()
    out = (tmp**2).mean()
    return out

t = round(dat[:,0].mean())
kb =  8.617e-5 # eV/K
kb_natoms = kb*args.num


cv_ion = fluct(dat[args.beg:,2])/(t**2)/kb/kb_natoms
tmp = []
for i in range(len(dat)):
    cv_i = fluct(dat[i:,2])/(t**2)/kb/kb_natoms
    tmp.append(cv_i)

s_ele = -dat[:,1].mean()/t
print('Temperature K |', 'Volume A3 |', 'Cv ionic Nkb |', 'Sele (eV/K)')
print('{0} \t {1} \t {2} \t {3}'.format(t,vol,cv_ion,s_ele))

fig,ax = plt.subplots(2,1,figsize=(6,6),sharex=True,sharey=False)
ax[0].plot(dat[:,2])
ax[1].plot(tmp)
ax[0].grid()
ax[0].minorticks_on()
ax[1].set_xlabel('step')
ax[0].set_ylabel('Eng without entropy')
ax[1].set_ylabel('cv_ion')
ax[1].minorticks_on()
ax[1].grid()
plt.show()



