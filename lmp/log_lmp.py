#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 00:21:59 2020
modified from cmd_interface.py

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
parser.add_argument("input_file", type=str, help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-x", type=str, default="Step", help="Data to plot on the first axis")
parser.add_argument("-y", type=str, nargs="+", help="Data to plot on the second axis. You can supply several names to get several plot lines in the same figure.")
parser.add_argument("-a", "--running_average", default=False, action='store_true', help="average y?")
parser.add_argument("-r", "--run_num", type=int, default=-1, help="run_num should be set if there are several runs and thermostyle does not change from run to run")
parser.add_argument("-s", "--store", default=False, action='store_true', help="Defualt:  Do not save data as outfile")
parser.add_argument("-of", "--outfile",type=str,default='log.properties', help="out file name")
parser.add_argument("-np", "--no_plot", default=False, action='store_true', help="Defualt:  Do not plot")


args = parser.parse_args()

log = File(args.input_file)
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


plt.figure()
for i in range(len(args.y)):
    plt.plot(x,ys[i],label=args.y[i])
if args.running_average:
    average=np.array(ys).mean(axis=0)
    plt.plot(x,average,label='average')
plt.legend()
plt.grid(True)
plt.show()


if args.store:
    header=args.x + ' '+  ' '.join(args.y)
    np.savetxt(args.outfile,np.concatenate(([x],ys),axis=0).T,fmt='%12.10f',header = header)
