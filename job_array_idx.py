#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 12:26:28 2020
https://www.geeksforgeeks.org/python-different-ways-to-kill-a-thread/
TODO: add condition to stop running thread after reach the last step => done
1-submit vasp job based on lammps dump file
2-stop at max job id prescribed
@author: jiedeng
"""

import threading
from subprocess import call
import os
import datetime
import argparse
import numpy as np
import sys
sys.setrecursionlimit(100000)
parser = argparse.ArgumentParser()
parser.add_argument("--max_job","-mj",type=int,help="max job allowed")
parser.add_argument("--idx","-id",type=str,help="index file, relative index only!")
#parser.add_argument("--count","-id",type=int,default=0,help="index file, relative index only!")

def get_waiting_jobs():
    call("/u/systems/UGE8.6.4/bin/lx-amd64/qstat -u jd848 | grep qw |wc|awk '{print $1}'>tmp.txt",shell=True)
    fp = open('tmp.txt')
    num_waiting_jobs = fp.readlines()
    fp.close()
    call("/u/systems/UGE8.6.4/bin/lx-amd64/qstat -u jd848 |wc|awk '{print $1}'>tmp.txt>tmp.txt",shell=True)
    fp = open('tmp.txt')
    num_tot_jobs = fp.readlines()
    fp.close()
    return int(num_waiting_jobs[0]), int(num_tot_jobs[0])-2

def get_max_job_id(path):
    recal=os.path.join(path,'recal')
    dirs=os.listdir(recal)
    dirs_int = [int(dd) for dd in dirs]
    max_job_id = 0
    for dd in dirs_int:
        if max_job_id<dd:
            max_job_id=dd
    return max_job_id,len(dirs_int)

## get the initial job index ##

tmp = 0 

args   = parser.parse_args()
idx    = np.loadtxt(args.idx).astype(int)

threshold  = 50  # waiting job threshold, if > threshold, do NOT submit, else, submit
next_batch = 20
#time_gap   = 100 # sec
job_lim    = 490 # in UCLA hoffmann
cwd = os.getcwd()
path = cwd
counti = 0
max_job_id,count0 = get_max_job_id(path)  
if max_job_id == idx[count0-1]:
    counti = count0  # start from counti

def foo(count=counti):
    
    num_waiting_jobs, num_tot_jobs =get_waiting_jobs()    
    now   = datetime.datetime.now()
    time  = now.strftime("%Y-%m-%d %H:%M:%S")
    time_gap = num_waiting_jobs # if num_waiting_jobs == 0, do not wait
     
    if num_waiting_jobs < threshold and (num_tot_jobs+next_batch) < job_lim:    
        max_job_id,_ = get_max_job_id(path)   
        print("{4}: # of waiting jobs: {0} < {1}, submit range({2},{3})".format(
              num_waiting_jobs,threshold, max_job_id, max_job_id+next_batch,time))
        for i in range(next_batch):
            nsw_range = '{0}-{1}'.format(idx[count],idx[count]+1)
#            print("python ~/script/mldp/recal_lmp.py -r {0}".format(nsw_range))
            call("python ~/script/mldp/recal_lmp.py -r {0}".format(nsw_range),shell=True)
            count = count+1
    elif num_waiting_jobs > threshold:
        print("{3}: # of waiting jobs: {0} > {1}, wait for {2} sec(s)".format(
              num_waiting_jobs,threshold, time_gap,time))
    elif (num_tot_jobs+next_batch) > job_lim:
        print("{3}: # of tot job {0} + # nextbatch > {1}, wait for {2} sec(s)".format(
              num_tot_jobs, job_lim, threshold, time_gap,time))
    max_job_id, _ = get_max_job_id(path)
    if  max_job_id < args.max_job:          
        threading.Timer(time_gap, foo(count)).start()
    else:
        quit()

foo()

