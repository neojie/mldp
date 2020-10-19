#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 09:11:02 2020

@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot(e):
    plt.figure()
    plt.plot(e[:,0],'o',label = 'data')
    plt.plot(e[:,1],'*',label = 'pred')
    plt.xlabel('step')
    plt.ylabel('TOTEN (eV)')
    if args.title:
        plt.title(args.title + ' atoms')
    plt.legend()
    plt.show()

e=np.loadtxt("result.e.out")

parser = argparse.ArgumentParser()
parser.add_argument("--title","-t",type = str,help="title")
parser.add_argument("--excudeupper","-eu",type = float,help="exclude data of energy > input value ")
parser.add_argument("--excudelower","-el",type = float,help="exclude data of energy < input value ")

args   = parser.parse_args()

if args.excudeupper:
    dat = e[e[:,0]<args.excudeupper]
if args.excudelower:
    try:
        dat = dat[dat[:,0]>args.excudelower]
    except:
        dat = e[e[:,0]>args.excudelower]
if args.excudeupper or args.excudelower:
    sqr_error = (dat[:,1] - dat[:,0])**2
    rmse = np.sqrt(np.mean(sqr_error))
    print("reivsed RMSE is", rmse)
    plot(dat)

else:
    plot(e)
    


