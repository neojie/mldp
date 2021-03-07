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


parser = argparse.ArgumentParser(description="Cv, assume run 0 is NVE")
parser.add_argument("--input_file","-i", type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-x", type=str, default="Step", help="Data to plot on the first axis")
parser.add_argument("-p", "--plot", default=False, action='store_true', help="Defualt: plot")
parser.add_argument("-n", "--natoms", default=160, type=int, help="natoms")
parser.add_argument("--volume", "-v",type=float,help='volume in A3')
parser.add_argument("--cv_nk", "-cv",type=float,help='cv in kb*N')


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
    
    
if args.volume:
    V = args.volume
else:
    infile = 'conf.lmp'
    fp = open(infile)
    ins = fp.readlines()
    print("  ?? vol not provided, parse from in conf.lmp") 
    xhi = False
    yhi = False
    zhi = False
    for line in ins:
        if 'xlo' in line:
            xlo = float(line.split()[0])
            xhi = float(line.split()[1])
        if 'ylo' in line:
            ylo = float(line.split()[0])
            yhi = float(line.split()[1])
        if 'zlo' in line:
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])
        if xhi and yhi and zhi:
            break
        
    V = (xhi - xlo)*(yhi - ylo)*(zhi - zlo)
    print(' ** volume = ', V)
    
# run_num  = 0
print('-----NVE-----')
run_num  = 0
x   = log.get(args.x,run_num=run_num)
args_y = ['Press','TotEng','Temp']
ys0  = [log.get(y,run_num=run_num) for y in args_y]


Step = log.get('Step',run_num=run_num)
if not check(Step):
#    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
    x = (x[selected_idx]).astype(float)
    ys0 = [(y[selected_idx]).astype(float) for y in ys0]
    print('**Fixed**')


kb =  8.617e-5 # eV/K
ys0 = np.array(ys0)
kb_natoms = kb*args.natoms

TotEng = ys0[1,:]
KinEng = 3/2*ys0[-1,:]*kb_natoms
PotEng = TotEng - KinEng
Temp   = ys0[-1,:]
Press  = ys0[-2,:]

xlen     = len(x)
xrange   = range(xlen//2,round(xlen/1.1)) # the last few alw
x        = x[xrange]
ys       = np.array(ys0)[:,xrange]
average0 =ys.mean(axis=1)
print_list(average0)

#tmp = ((ys - np.reshape(average0,(len(average0),1)))**2).mean(axis=1)
print("TotEng",xrange)
print(np.mean(TotEng))

print("KinEng",xrange)
print(np.mean(KinEng))

print("PotEng",xrange)
print(np.mean(PotEng))

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

pot   = PotEng[xrange]
temp  = average0[2]
out = []        
for i in range(len(pot)):
    cv_i = 1.5/(1-fluct(pot[i:])/1.5/kb/kb/args.natoms/temp/temp)
#    cv_i = fluct(pot[i:])/(temp**2)/kb/kb_natoms
    out.append(cv_i)

cv_mean = np.mean(out)
print('cv mean is',cv_mean)


deltaP = Press[xrange] - np.mean(Press[xrange])
deltaPot = PotEng[xrange] - np.mean(PotEng[xrange])


j2eV = 6.242e18
if args.cv_nk:
    term2 = 2/3*args.cv_nk*args.natoms*kb/V
else:  
    term2 = 2/3*cv_mean*args.natoms*kb/V
ev_A3_to_Pa = 1/j2eV/(1e-30)
gamma = [] # dP/dT
for i in range(len(deltaP)):
    
    term1 =(1-(deltaP[:(i+1)]*1e5*deltaPot[:(i+1)]).mean()*V*1e-30*j2eV/args.natoms/kb/kb/temp/temp)    
    gamma.append(term1*term2*ev_A3_to_Pa)

fig,ax = plt.subplots(3,1,figsize=(6,6),sharex=True,sharey=False)
ax[0].plot(pot)
ax[1].plot(out)
ax[2].plot(gamma)
ax[0].grid()
ax[1].grid()
ax[2].grid()
ax[0].minorticks_on()
ax[1].set_xlabel('step')
ax[0].set_ylabel('PotEng')
ax[1].set_ylabel('cv_ion')
ax[1].minorticks_on()
plt.show()


#plt.plot(ys[1] - average0[1],'*-');plt.xlabel('step');plt.ylabel('TotEng - mean');plt.show()
#plt.plot(ys[-2] - average0[-2],'*-');plt.xlabel('step');plt.ylabel('PotEng - mean');plt.show()
#plt.plot(ys[-1] - average0[-1],'*-');plt.xlabel('step');plt.ylabel('KinEng - mean');plt.show()

