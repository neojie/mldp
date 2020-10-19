#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:14:45 2019

1- allow running in the local folder
2- if exist post_recal_done flag, skip
3- if not, check

@author: jiedeng
"""

from shared_functions import check_outcar_done, check_outcar_done_slow, load_paths
import argparse
import os
import glob
from subprocess import call

print("=="*40)
print(" "*20, "MUST HAVE DSQ MODULE AVAIL!!!", " "*20)
print(" "*20, "Make sure running job not overlapped with here", " "*20)
print("=="*40)

print()
call('squeue -u jd848 --format %Z > running', shell=True)
fp=open('running');running_jobs = fp.readlines();fp.close()
#print(running_jobs)
os.remove('running')

existing_jobs = {}
for job in running_jobs:
    if 'recal' in job:
        existing_job = job.replace('\n','')
        existing_jobs[existing_job] = None
    
print("exist jobs are", existing_jobs)

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    tmp = load_paths(inputpath)
    paths = [os.path.join(path,'recal') for path in tmp]
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd)
    paths = [cwd]

def build_dsq_file(path):

    done   = 0
    fail   = 0
    out    = []    

    subfolders = [f.path for f in os.scandir(path) if f.is_dir() ]
    
    for folder in subfolders:
        outcar = os.path.join(folder,'OUTCAR')
        flag = False
        if os.path.exists(outcar):
            try:
                outcar_done = check_outcar_done(outcar)
            except:
                outcar_done = check_outcar_done_slow(outcar)
            if outcar_done:
                flag = True               
        if flag:
            done += 1
        else:
            fail += 1
            target_folder = cwd
            failpath = os.path.join(target_folder,folder)     
            if fail < 1000:
                out.append('cd '+failpath+'; ' + 'mpirun vasp_std' +'\n')
            else:
                raise ValueError("job cannot exceed 1000!")

    fw   = open(os.path.join(path,'out'),'w')
    fw.writelines(out)
    fw.close()
    return done, fail
    
  
def remove_file(path,wild_card):
    target =  os.path.join(path,wild_card)
    for file_to_remove in glob.glob(target):
        os.remove(file_to_remove)

    
for path in paths:
    print(path)
    remove_file(path, '*.out')
    if os.path.exists(os.path.join(path,'post_recal_done')) or path in existing_jobs:
        print('### past recal done')
    else:
        done, fail = build_dsq_file(path)
        print("----"*10)
        print("There are {0} runs done and {1} fail".format(done, fail))
        print("----"*10)
        if fail == 0:
            print('### all past recal done')
            with open(os.path.join(path,'post_recal_done'), 'w') as fp: 
                pass
            print('### remove all useless files')  ## to do later
            remove_file(path, 'out*');remove_file(path, '*sh');remove_file(path, '*tsv')
        else:
            os.chdir(path)
            call("dsq --job-file out --cpus-per-task 48 --mem-per-cpu 6g -t 02:00:00 -p mpi --mail-type None --batch-file sub_script.sh\n",shell=True)
            call("sbatch sub_script.sh", shell=True)            
            os.chdir(cwd)

