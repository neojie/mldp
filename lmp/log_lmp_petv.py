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


parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
parser.add_argument("--input_file","-i", type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-x", type=str, default="Step", help="Data to plot on the first axis")
parser.add_argument("-p", "--plot", default=True, action='store_false', help="Defualt: plot")


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
args_y = ['Press','TotEng','Temp', 'Volume']
ys  = [log.get(y,run_num=run_num) for y in args_y]

Step = log.get('Step',run_num=run_num)
if not check(Step):
#    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
    x = (x[selected_idx]).astype(float)
    ys = [(y[selected_idx]).astype(float) for y in ys]
    print('**Fixed**')

average0=np.array(ys).mean(axis=1)
print_list(average0)
if args.plot:   
    fig,ax = plt.subplots(3,2,figsize=(10,6),sharex=False,sharey=False)
    for i in range(3):
        ax[i][0].plot(x,ys[i],label=args_y[i])
        ax[i][0].legend()
        ax[i][0].grid(True)


print('-----NVE-----')
run_num  = 1
x   = log.get(args.x,run_num=run_num)
args_y = ['Press','TotEng','Temp']
ys  = [log.get(y,run_num=run_num) for y in args_y]

Step = log.get('Step',run_num=run_num)
if not check(Step):
#    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
    x = (x[selected_idx]).astype(float)
    ys = [(y[selected_idx]).astype(float) for y in ys]
    print('**Fixed**')

average1=np.array(ys).mean(axis=1)
print_list(average1)

print('-----NVT + NVE-----')
print_list(np.concatenate(average0 , average1))

if args.plot:
    for i in range(3):
        ax[i][1].plot(x,ys[i],label=args_y[i])
        ax[i][1].legend()
        ax[i][1].grid(True)

    plt.show()
