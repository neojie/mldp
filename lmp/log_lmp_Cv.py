#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 00:21:59 2020
modified from cmd_interface.py

WARNING message can only be handled if filed is less than lammps fields
say if I have
Step TotEng
then this code cannot handle it. => so I use 'sed' as workaround



----  modification of log_lmp.py for output specifically
1- NVT run : pressure,   TotEng, temperature, volume
2- NVE run : pressure, temperature, , TotEg
@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
#import os
import argparse
#import glob
#from util import blockAverage
from lammps_logfile import File


parser = argparse.ArgumentParser(description="Cv, assume run 0 is NVT")
parser.add_argument("--input_file","-i", type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-x", type=str, default="Step", help="Data to plot on the first axis")
parser.add_argument("-p", "--plot", default=False, action='store_true', help="Defualt: plot")
parser.add_argument("-n", "--natoms", default=160, type=int, help="natoms")


args = parser.parse_args()

try:
    log = File(args.input_file)
except:
    from subprocess import call
    call("cp {0} {1}".format(args.input_file, args.input_file+'tmp'), shell=True) # do not change in the original file, better for checking on running sinulation
    call("sed -i 's/style restartinfo set but has//' {0}".format(args.input_file+'tmp'), shell=True)
    log = File(args.input_file+'tmp')
    call("rm {0}".format(args.input_file+'tmp'), shell=True)


"""
log file may have WARNING message mess everything up
>>> log.get("Temp")
array(['0', 'Pair', '0', ..., '0', 'Pair', '0'], dtype=object)
Pyz=log.get('Pyz')
array([ 13064.555  ,          nan,   -269.48088,   -269.48088])
"""
def check(dat):
    for ele in dat:
        if str(ele).replace('.','').replace('e','').replace('d','').replace('+','').replace('-','').isdigit(): # data may be in int or float, may contain . , e , d , + , -
            pass
        else:
            return False
    return True

def select(dat):
    selected_idx = []
    for i  in range(len(dat)):
        if str(dat[i]).isdigit():
            selected_idx.append(i)
    return selected_idx

def print_list(L):
    print('\t'.join([str(i) for i in L]))
# run_num  = 0
print('-----NVT-----')
run_num  = 0
x   = log.get(args.x,run_num=run_num)
args_y = ['Press','TotEng','Temp', 'PotEng','KinEng']
ys0  = [log.get(y,run_num=run_num) for y in args_y]


Step = log.get('Step',run_num=run_num)
if not check(Step):
#    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
    x = (x[selected_idx]).astype(float)
    ys0 = [(y[selected_idx]).astype(float) for y in ys0]
    print('**Fixed**')

xlen = len(x)
xrange = range(xlen//2,xlen)
x  = x[xrange]
ys = np.array(ys0)[:,xrange]
average0=ys.mean(axis=1)
print_list(average0)

kb =  8.617e-5 # eV/K
kb_natoms = kb*args.natoms

tmp = ((ys - np.reshape(average0,(len(average0),1)))**2).mean(axis=1)
print("TotEng",xrange)
print(tmp[1]/(average0[2]**2)/kb/kb_natoms)

print("KinEng",xrange)
print(tmp[-1]/(average0[2]**2)/kb/kb_natoms)

print("PotEng",xrange)
print(tmp[-2]/(average0[2]**2)/kb/kb_natoms)

if args.plot:   
    fig,ax = plt.subplots(3,2,figsize=(10,6),sharex=False,sharey=False)
    for i in range(3):
        ax[i][0].plot(x,ys[i],label=args_y[i])
        ax[i][0].legend()
        ax[i][0].grid(True)

def fluct(data):
    tmp=data - data.mean()
    out = (tmp**2).mean()
    return out

pot = np.array(ys0)[-2,:]
out = []        
for i in range(len(pot)):
    cv_i = fluct(pot[i:])/(average0[2]**2)/kb/kb_natoms
    out.append(cv_i)

fig,ax = plt.subplots(2,1,figsize=(6,6),sharex=True,sharey=False)
ax[0].plot(pot)
ax[1].plot(out)
ax[0].grid()
ax[0].minorticks_on()
ax[1].set_xlabel('step')
ax[0].set_ylabel('PotEng')
ax[1].set_ylabel('cv_ion')
ax[1].minorticks_on()
ax[1].grid()
plt.show()


#plt.plot(ys[1] - average0[1],'*-');plt.xlabel('step');plt.ylabel('TotEng - mean');plt.show()
#plt.plot(ys[-2] - average0[-2],'*-');plt.xlabel('step');plt.ylabel('PotEng - mean');plt.show()
#plt.plot(ys[-1] - average0[-1],'*-');plt.xlabel('step');plt.ylabel('KinEng - mean');plt.show()

