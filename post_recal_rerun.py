#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 18:48:47 2020

*** to rerun sims that didn't start or some issue with OUTCAR (didn't finish, converge ...) ***
Updated version of {check_nbands_nelm.py + post_recal.py}
Example command: python $mldp/post_recal_rerun.py -ip all -v -ss address/to/sub_vasp.sh > test

@author: akashgpt
"""

from shared_functions import load_paths, check_outcar_done, check_outcar_done_slow
import os
import argparse
cwd    = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file.\
                    1-None, current folder\
                    2-all, all subfolders \
                    3-file, folders stored at file")

parser.add_argument("--OUTCAR","-o",type = str, default = 'OUTCAR', help="OUTCAR name, if not 'OUTCAR', please specify the name")
#parser.add_argument("--verbose","-v",type = bool, default = True, help="OUTCAR name")

parser.add_argument('--verbose',"-v", default=True, action='store_false')
parser.add_argument('--remove_outcar',"-ro", default=False, action='store_true',help="removing all OUTCAR files that are not good")
parser.add_argument("--submissionscript","-ss",help="input path file")

args   = parser.parse_args()

if not args.inputpath:
    print("No folders point are provided. Use current working path")
    inputpath = os.path.join(cwd)
    paths = [cwd]
else:
    if args.inputpath == 'all':
        paths = []
        print("check all sub folders at {0}".format(cwd))
        for directory in os.listdir(cwd):
            if os.path.isdir(directory):
                paths.append(directory)
    else:
        print("Check files in {0}  ".format(args.inputpath))
        inputpath = args.inputpath
        paths = load_paths(inputpath)        


def check(path):
    outcar = os.path.join(path,args.OUTCAR)
    fp=open(outcar)
    flag_nelm, flag_nbands, flag_iter, flag_occ = False, False, False, False
    check_iter, check_nbands = False, False
    good_iter, good_nbands = True, False
    outcar_done_flag = False

    if os.path.exists(outcar):
        try:
            outcar_done = check_outcar_done(outcar)
        except:
            outcar_done = check_outcar_done_slow(outcar)
        if outcar_done:
            outcar_done_flag = True
        # else:
        #     print("outcar not done",outcar)

    for i in range(1000000):
        line = fp.readline()
        if 'NELM' in line:#1st appear
            nelm = int(line.split()[2].replace(';','').replace(',','').replace('=',''))
            flag_nelm = True
        if 'number of bands    NBANDS' in line:#2nd appear
            nbands = int(line.split()[-1].replace(';','').replace(',','').replace('=',''))
            flag_nbands = True
        if 'Iteration' in line:#3rd appear
            iteration = int(line.split('(')[-1].split(')')[0])
            flag_iter = True
        if 'band No.' in line:
            for _ in range(nbands):
                line_kpoints = fp.readline()
#            print(line_kpoints)
            occupancy=float(line_kpoints.split()[-1])
            flag_occ = True
        if flag_nelm and flag_iter: 
            check_iter = True
            if nelm == iteration: ## iteration varies from [1 to nelm]
                good_iter = False
                
        if flag_nbands and flag_occ:
            check_nbands = True
            if occupancy<1e-6:
                good_nbands = True
                
        if check_iter and check_nbands:
            if args.verbose:
                print("details")
                print(nelm, iteration)
                print(nbands, occupancy)
            return good_iter, good_nbands, outcar_done_flag

    return good_iter, good_nbands, outcar_done_flag



rerun    = []  

print("\n\nKey:\n l0: OUTCAR does not exist;\n l1: OUTCAR exists BUT simulation never finished;\n l2: OUTCAR exists and simulation finished BUT nbands not enough;\n l3: OUTCAR exists, simulation finished and nbands enough BUT # of iterations not enough for convergence\n\n")
print("[Folder Path] [Flag]")

for path in paths:

    outcar = os.path.join(path,args.OUTCAR)

    if not os.path.exists(outcar): # OUTCAR does not exist
        print(path, "l0")#if OUTCAR does not exist
        target_folder = cwd
        failpath = os.path.join(target_folder,path)
        rerun.append('cd '+failpath+'; ' + 'sbatch {0}'.format(args.submissionscript) +'\n')

    if os.path.exists(outcar): # OUTCAR exists
        good_iter, good_nbands, outcar_done_flag  = check(path)
        # if not good_iter:
        #     print("bad iter",path)

        if not outcar_done_flag: # OUTCAR not done/simulation never finished
            print(path, "l1")#if OUTCAR exists BUT simulation never finished
            target_folder = cwd
            failpath = os.path.join(target_folder,path)
            rerun.append('cd '+failpath+'; ' + 'sbatch {0}'.format(args.submissionscript) +'\n')

        if outcar_done_flag: # if OUTCAR completed/simulation finished
            if not good_nbands: # if nbands not enough
                print(path, "l2")#if OUTCAR exists and simulation finished BUT nbands not enough
                # print("bad nbands",path)
                target_folder = cwd
                failpath = os.path.join(target_folder,path)
                rerun.append('cd '+failpath+'; ' + 'sbatch {0}'.format(args.submissionscript) +'\n')

            if good_nbands: # if nbands enough
                if not good_iter: # if number of iterations not enough for convergence
                    print(path, "l3")#if OUTCAR exists, simulation finished and nbands enough BUT iteration not enough
                    target_folder = cwd
                    failpath = os.path.join(target_folder,path)
                    rerun.append('cd '+failpath+'; ' + 'sbatch {0}'.format(args.submissionscript) +'\n')




        if args.remove_outcar and ((not good_iter) or (not good_nbands)):
            os.remove(os.path.join(path,args.OUTCAR))
            print("** remove outcar successfully")
    else:
        if args.verbose:
            print("NO OUTCAR")

# print()
fw   = open(os.path.join(cwd,'rerun'),'w')
rerun.append('cd ..')
fw.writelines(rerun)
fw.close()

print("\n\nRun 'source rerun' to rerun the jobs that did not finish properly.")
            
