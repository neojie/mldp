#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 12:26:28 2020
https://www.geeksforgeeks.org/python-different-ways-to-kill-a-thread/
TODO: add condition to stop running thread after reach the last step => done
1-submit vasp job based on out file generated by post_recal_lmp.py
2-out file can be run by `bash out`
3-maintain that total number of job waiting and running < 490
@author: jiedeng
"""
import sys
sys.setrecursionlimit(10000)
import threading
from subprocess import call
import os
import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--out","-o",type=str,help="out file path")
parser.add_argument("--outs","-os",type=str,help="file of out files' paths")

args   = parser.parse_args()


if args.out:
    fp=open(args.out)
    jobs = fp.readlines()
if args.outs:
    fp=open(args.outs)
    jobs = []
    for out in fp:
        tmp = open(out.split()[0])
        for tmp2 in tmp.readlines():
            jobs.append(tmp2)
print('There are {0} jobs.'.format(len(jobs)))
print('The first jobs is', jobs[0])
print('The last jobs is', jobs[-1])

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

def foo(idx):
    num_waiting_jobs, num_tot_jobs =get_waiting_jobs()    
    now   = datetime.datetime.now()
    time  = now.strftime("%Y-%m-%d %H:%M:%S")
    time_gap = num_waiting_jobs # if num_waiting_jobs == 0, do not wait
    
    if num_waiting_jobs < threshold and (num_tot_jobs+next_batch) < job_lim:    
        print("{2}: # of waiting jobs: {0} < {1})".format(
              num_waiting_jobs,threshold,time))
    
        print("The first job is", jobs[idx:idx+next_batch][0])
        print("The last job is", jobs[idx:idx+next_batch][-1])
        for job in jobs[idx:idx+next_batch]:
            call(job,shell=True)
        idx = idx+next_batch
    elif num_waiting_jobs > threshold:
        print("{3}: # of waiting jobs: {0} > {1}, wait for {2} sec(s)".format(
              num_waiting_jobs,threshold, time_gap,time))
    elif (num_tot_jobs+next_batch) > job_lim:
        print("{3}: # of tot job {0} + # nextbatch > {1}, wait for {2} sec(s)".format(
              num_tot_jobs, job_lim, threshold, time_gap,100))
    if  idx < len(jobs):          
        threading.Timer(time_gap, foo(idx)).start()
    else:
        quit()

threshold  = 50  # waiting job threshold, if > threshold, do NOT submit, else, submit
next_batch = 10
job_lim    = 490 # in UCLA hoffmann
cwd = os.getcwd()
path = cwd
foo(0)

