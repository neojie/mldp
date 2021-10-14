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
import sys
import time
sys.setrecursionlimit(1000000)
recal_path = os.path.join(os.getcwd(),'recal')
try:
    os.mkdir(recal_path)     
except:
    print('***recal exists in',recal_path)
    
parser = argparse.ArgumentParser()
parser.add_argument("--max_job","-mj",type=int,help="max job allowed, typically set to be number of frames in the deepmd")
parser.add_argument("--deepmd","-d",help="deepmd path file")
parser.add_argument("--inputfile","-if",help="input files for vasp cal, default is cwd+inputs, if input, please input the absolute path")
args   = parser.parse_args()
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
    return max_job_id

def foo():
    num_waiting_jobs, num_tot_jobs =get_waiting_jobs()    
    now   = datetime.datetime.now()
    timex  = now.strftime("%Y-%m-%d %H:%M:%S")
#    time_gap = num_waiting_jobs # if num_waiting_jobs == 0, do not wait
     
    if num_waiting_jobs < threshold and (num_tot_jobs+next_batch) < job_lim:    
        max_job_id = get_max_job_id(path)
        if max_job_id+next_batch>args.max_job:
            end_job = args.max_job+1
        else:
            end_job = max_job_id+next_batch
        print("{4}: # of waiting jobs: {0} < {1}, submit range({2},{3})".format(
              num_waiting_jobs,threshold, max_job_id, max_job_id+next_batch,timex))
        nsw_range = '{0}-{1}'.format(max_job_id,end_job)
        call("python ~/script/mldp/recal_dpdata.py -r {0} -d {1} -if {2}".format(
                nsw_range,args.deepmd, args.inputfile),shell=True)
    elif num_waiting_jobs > threshold:
        print("{3}: # of waiting jobs: {0} > {1}, wait for {2} sec(s)".format(
              num_waiting_jobs,threshold, num_waiting_jobs, timex))
        time.sleep(num_waiting_jobs)
    elif (num_tot_jobs+next_batch) > job_lim:
        print("{3}: # of tot job {0} + # nextbatch > {1}, wait for {2} sec(s)".format(
              num_tot_jobs, job_lim, threshold, num_waiting_jobs,timex))
        time.sleep(threshold)
    max_job_id = get_max_job_id(path)
    if  max_job_id < args.max_job:          
        threading.Timer(num_waiting_jobs, foo).start()
    else:
        quit()



threshold  = 450  # waiting job threshold, if > threshold, do NOT submit, else, submit
next_batch = 20
#time_gap   = 100 # sec
job_lim    = 480 # in UCLA hoffmann, actually 500, intentionally leaving 20 
cwd = os.getcwd()
path = cwd
foo()

