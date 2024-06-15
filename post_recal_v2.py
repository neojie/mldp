#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:14:45 2019, Modified on 20240601

1- allow running in the local folder
2- if exist post_recal_done flag, skip
3- if not, check

*** to check if all runs done successfully ***
Now modified to work on Princeton clusters
Example command: python $mldp/post_recal_v2.py > checks && python $mldp/check_nbands_nelm.py -ip all -v >> check

@author: jiedeng, akashgpt
"""

from shared_functions import check_outcar_done, check_outcar_done_slow
import argparse
import os
import glob
from subprocess import call

print("=="*40)
#print(" "*20, "MUST HAVE DSQ MODULE AVAIL!!!", " "*20)
print(" "*20, "Running jobs are not considered", " "*20)
print("=="*40)

print()
call("squeue -u $USER | grep qw |wc|awk '{print $1}'>tmp.txt",shell=True)
fp=open('tmp.txt');running_jobs = fp.readlines();fp.close()
# print(running_jobs)
# os.remove('running')

existing_jobs = {}
for job in running_jobs:
    if 'recal' in job:
        existing_job = job.replace('\n','')
        existing_jobs[existing_job] = None
    
print("exist jobs are", existing_jobs)

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-if",help="where the recal folder is, default is cwd")
parser.add_argument("--submissionscript","-ss",help="input path file")

args   = parser.parse_args()

cwd    = os.getcwd()

if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
#    inputpath = args.inputpath
#    tmp = load_paths(inputpath)
    paths = [args.inputpath]
else:
    print("No folders are provided. Use default value folders")
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
            # else:
            #     print(folder, "flag false")               
        if flag:
            done += 1
        else:
            fail += 1
            print(folder)
            target_folder = cwd
            failpath = os.path.join(target_folder,folder)     
            out.append('cd '+failpath+'; ' + 'sbatch {0}'.format(args.submissionscript) +'\n')

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
            print("To re-submit the failed jobs use post_recal_rerun.py")
#            os.chdir(path)
#            call("dsq --job-file out -c 4 --mem-per-cpu 6g -t 08:00:00 -p scavenge --mail-type None --batch-file sub_script.sh\n",shell=True)
#            call("sbatch sub_script.sh", shell=True)            
#            os.chdir(cwd)

