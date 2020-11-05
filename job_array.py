#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 12:26:28 2020
https://www.geeksforgeeks.org/python-different-ways-to-kill-a-thread/

TODO: add condition to stop running thread after reach the last step
@author: jiedeng
"""

import threading
from subprocess import call
import os
import datetime

def get_waiting_jobs():
    call("qstat -u jd848 | grep qw |wc|awk '{print $1}'>tmp.txt",shell=True)
    fp = open('tmp.txt')
    num_waiting_jobs = fp.readlines()
    fp.close()
    call("qstat -u jd848 |wc|awk '{print $1}'>tmp.txt>tmp.txt",shell=True)
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
    time  = now.strftime("%Y-%m-%d %H:%M:%S")
    if num_waiting_jobs < threshold and (num_tot_jobs+next_batch) < job_lim:    
        max_job_id = get_max_job_id(path)
        
        print("{4}: # of waiting jobs: {0} < {1}, submit range({2},{3})".format(
              num_waiting_jobs,threshold, max_job_id, max_job_id+next_batch,time))
        nsw_range = '{0}-{1}'.format(max_job_id,max_job_id+next_batch)
        call("python ~/script/mldp/recal_lmp.py -r {0}".format(nsw_range),shell=True)
    elif num_waiting_jobs < threshold:
        print("{3}: # of waiting jobs: {0} > {1}, wait for another {2} sec(s)".format(
              num_waiting_jobs,threshold, time_gap,time))        
    threading.Timer(time_gap, foo).start()

threshold = 50
next_batch = 20
time_gap = 100
job_lim = 500
cwd = os.getcwd()
path = cwd
foo()

