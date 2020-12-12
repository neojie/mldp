#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:38:32 2020

@author: jiedeng
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nsw","-n",type=int,help="nsw of OUTCAR")
parser.add_argument("--n_idx","-ni",type=int,help="number of idx")
args   = parser.parse_args()

x = np.random.choice(range(args.nsw), args.n_idx, replace=False)
x.sort()
np.savetxt('id_rand_{0}'.format(args.n_idx),x,fmt='%d')