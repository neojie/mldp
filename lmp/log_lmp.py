#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 00:21:59 2020
modified from cmd_interface.py

WARNING message can only be handled if filed is less than lammps fields
say if I have
Step TotEng
then this code cannot handle it. => so I use 'sed' as workaround

@author: jiedeng
"""


import argparse

parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
#parser.add_argument("input_file", type=str, help="Lammps log file containing thermo output from lammps simulation.") # default cannot be set to positional argument
parser.add_argument("--input_file",'-i', type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-x", type=str, default="Step", help="Data to plot on the first axis")
parser.add_argument("-y", type=str, nargs="+", help="Data to plot on the second axis. You can supply several names to get several plot lines in the same figure.")
parser.add_argument("-a", "--running_average", default=True, action='store_false', help="Default: average y and print out the averaged value ")
parser.add_argument("-h", "--half_window", default=True, action='store_false', help="Default: average the second half of y")
parser.add_argument("-r", "--run_num", type=int, default=-1, help="run_num should be set if there are several runs and thermostyle does not change from run to run")
parser.add_argument("-s", "--store", default=False, action='store_true', help="Defualt:  Do not save data as outfile")
parser.add_argument("-of", "--outfile",type=str,default='log.properties', help="out file name")
parser.add_argument("-p", "--plot", default=True, action='store_false', help="Defualt: plot")


args = parser.parse_args()

from lammps_logfile import File
import numpy as np
import matplotlib.pyplot as plt

try:
    log = File(args.input_file)
except:
    from subprocess import call
    call("cp {0} {1}".format(args.input_file, args.input_file+'tmp'), shell=True) # do not change in the original file, better for checking on running sinulation
    call("sed -i 's/style restartinfo set but has//' {0}".format(args.input_file+'tmp'), shell=True)
    log = File(args.input_file+'tmp')
    call("rm {0}".format(args.input_file+'tmp'), shell=True)
    
x   = log.get(args.x,run_num=args.run_num)
ys  = [log.get(y,run_num=args.run_num) for y in args.y]
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


Step = log.get('Step',run_num=args.run_num)
if not check(Step):
    print('**data messed up**')
    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
    x = (x[selected_idx]).astype(float)
    ys = [(y[selected_idx]).astype(float) for y in ys]
    print('**Fixed**')
    
if args.running_average:
    if args.half_window:
        xrange = range(len(x)//2,len(x))
    else:
        xrange = range(len(x))
    print('----average of step {0}----'.format(xrange))
    print(args.y)
    if len(args.y) == 1:
        ys_mean = [ys[0][xrange].mean()]
    else:
        ys_mean = [y[xrange].mean(axis=0) for y in ys]
    print('\t'.join([str(i) for i in ys_mean]))
    
if args.plot:
    horizontal_scale = len(args.y)*3
    fig,ax = plt.subplots(len(args.y),1,figsize=(6,horizontal_scale),sharex=True,sharey=False)
    if len(args.y) == 1:
        ax.plot(x,ys[0],label=args.y[0])
        ax.legend();ax.grid()
    else:
        for i in range(len(args.y)):
            ax[i].plot(x,ys[i],label=args.y[i])
            ax[i].legend();ax[i].grid()
    plt.show()


if args.store:
    header=args.x + ' '+  ' '.join(args.y)
    np.savetxt(args.outfile,np.concatenate(([x],ys),axis=0).T,fmt='%12.10f',header = header)
