#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:02:21 2019
version 2 
#readlines behave very funny, cannot understand...
# first build folder with XDATCAR, INCAR,POTCAR, bs, KPOINTS


@author: jiedeng
"""

import numpy as np
import os
from shutil import copy
import re
from shared_functions import load_paths, reverse_readline, check_outcar_done, check_incar

def xdt2pos(path,sel_nsw,tot_ele,copybs=False):
    """
    build POSCAR based on XDATCAR given in the path
    copy INCAR, POTCAR, KPOINTS into the same folder
    
    """

    XDATCAR = open(os.path.join(path,"XDATCAR"),'r')
    incar = os.path.join(path,'INCAR')
    path = os.path.join(path,'recal')

    for _ in range(7):
        XDATCAR.readline()
    assert len(sel_nsw)>1  ## at least two iteration are selected

    skp_nsw  = []
    for i in range(len(sel_nsw)):
        if i == 0:
            skp_nsw.append(sel_nsw[i]*(tot_ele+1))
        else:
            skp_nsw.append((sel_nsw[i] - sel_nsw[i-1] -1)*(tot_ele+1))
    for i,diff_i in zip(sel_nsw,skp_nsw):
        try:
            os.mkdir(os.path.join(path,str(i+1))) # return None
        except:
            print("Folder {0} already exists,skip making".format(i))
        target_path = os.path.join(path,str(i+1))
#        print(target_path)
        ls_target_path = os.listdir(target_path)
                # skip some portion
        for _ in range(diff_i):
            XDATCAR.readline()
            
        if 'OUTCAR'in ls_target_path and check_outcar_done(os.path.join(target_path,'OUTCAR')):
            print('calculation done')
            for _ in range(tot_ele+1):       
                XDATCAR.readline()
                
        elif 'INCAR'   in ls_target_path and \
             'POTCAR'  in ls_target_path and \
             'POSCAR'  in ls_target_path and \
             'KPOINTS' in ls_target_path:
            for _ in range(tot_ele+1):       
                XDATCAR.readline()
                
        else:                
            coord = XDATCAR.readline() # Direct
#            print(coord)
            atomic_position = []
            for j in range(tot_ele):       
                line_tmp = XDATCAR.readline()
                atomic_position.append(line_tmp)
            
            tobewriten = []
            tobewriten.append(title+'\n')
            tobewriten.append('{0}\n'.format(scaling_factor))
            for j in range(3):
                tobewriten.append('    {0:14.9f}{1:14.9f}{2:14.9f}\n' .format(lattice[0][j],lattice[1][j],lattice[2][j]))
            for j in range(len(list_ele)):
                tobewriten.append('    %4s%s'%(list_ele[j],' '))
            tobewriten.append('\n')
            for j in range(len(num_ele)):
                tobewriten.append('    %4d%s'%(num_ele[j],' '))   
            tobewriten.append('\n')
            tobewriten.append('Direct\n')
            for j in range(tot_ele):
                tobewriten.append(atomic_position[j])    
            #begin to write POSCAR content to tmp file
            fw = open(os.path.join(target_path,'POSCAR'),'w')
            fw.writelines(tobewriten)
            fw.close()  
#            os.rename('tmp',os.path.join(target_path,'POSCAR'))
            target_incar = check_incar(incar)
            copy(os.path.join(inputfile, target_incar),target_path)
            os.rename(os.path.join(target_path,target_incar),os.path.join(target_path,'INCAR'))
            copy(os.path.join(inputfile, 'KPOINTS'),target_path)
            copy(os.path.join(inputfile, 'POTCAR'),target_path)

            
#            if copybs:
#                copy('bs.sh',target_path)
#                copy('bs_job.sh',target_path)
    XDATCAR.close()

############################ read from file for paths  ############################
### get target path


def check_xdatcar_last_frame(xdatcar,tot_ele):
    """
    extract the last frame
    extract tot NSW
    """
    last_frame = []
    XDATCAR_reverse = reverse_readline(xdatcar)
    for _ in range(tot_ele):
        line =  next(XDATCAR_reverse)
        last_frame.append([float(j) for j in line.split()])
    tot_nsw  = int(re.sub("[^0-9]", "",next(XDATCAR_reverse)))
    last_frame.reverse()
    return tot_nsw, np.array(last_frame)



def create_job(path,tot_nsw,tot_ele,is_liquid,liquid_frames = 500, solid_frames = 100,copybs = True):
        
    """
    output sel_nsw MUST be sorted!
    """
    # 1- get NSW
#    nsw = build_POSCARs(path)  #
    # 2- decide NSW to feed in
    #  2.1 judge if solid/liquid/glass
    #      solid = atoms move < 30% accross the whole iteration
    #      liquid and glass = atoms move > 30% accross the whole iteration
    #  2.2 liquid = 500 frame
    #      solid  = 100 frame
    
    if is_liquid:
#        sel_nsw = list(range(int(tot_nsw*0.05),tot_nsw,tot_nsw//liquid_frames)) # skip the first few percent
        sel_nsw = np.random.choice(range(tot_nsw), liquid_frames, replace=False)
    else:
        sel_nsw = np.random.choice(range(tot_nsw), solid_frames, replace=False)
    sel_nsw.sort()
#        sel_nsw = list(range(int(tot_nsw*0.05),tot_nsw,tot_nsw//solid_frames))
#    sel_nsw = [i+1 for i in sel_nsw]  # XDATCAR starting from 1, python starting from 0
    xdt2pos(path,sel_nsw,tot_ele,copybs = copybs)
    return sel_nsw
    

def atoms_moved(first_frame,last_frame,lattice,cutoff=0.3):
    
    diff  = first_frame - last_frame
    judge = np.logical_and(abs(diff)>cutoff, abs(diff)< (1-cutoff))
    if judge.any():
        print("is liquid")
        return True
    else:
        print("is solid")
        return False

def dsq_jobs(path):
    """
    to keep the set-up consistent, we still call this function as is
    But no dsq exist at ucla
    This function is used to sub_vasp.sh file only
    """
    target_folder = os.path.join(path,'recal')
    copy(os.path.join(inputfile, 'sub_vasp.sh'),target_folder)
      
        

from subprocess import call
def run_dsq(cwd,path,sel_nsw):
    target_folder = os.path.join(path,'recal')
    sub_vasp = os.path.join(target_folder,'sub_vasp.sh')
    for i in sel_nsw:
        jobpath = os.path.join(target_folder,str(i+1))  
        os.chdir(jobpath)
        call("sbatch {0}".format(sub_vasp), shell=True)
        os.chdir(cwd)
        

    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--inputfile","-if",help="input files")
parser.add_argument("--liquid_gap","-lg",type=int,default=1,help="data gap for liquid")
parser.add_argument("--solid_gap","-sg",type=int,default=1,help="data gap for solid")

args   = parser.parse_args()

#print(args.run_vasp)
cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
    paths = load_paths(inputpath)
else:
    print("input path are provided.")
    paths = [cwd]

if args.inputfile:
    print("Check files in {0}  ".format(args.inputfile))
    inputfile = args.inputfile
else:
    print("No folders point are provided. Use default value folders")
    inputfile = os.path.join(cwd,'inputs')
    
#paths = load_paths(inputpath)
for path in paths:
    print('###',path)
    # get basic info

    if os.path.exists(os.path.join(path,'recal')):
#        path = os.path.join(path,'recal')
#        pass 
        print('*** SKIP',path,'has recal folder')
    else:
        print('--> Build recal',path)
        os.mkdir(os.path.join(path,'recal'))        
        XDATCAR = open(os.path.join(path,"XDATCAR"),'r')
        title   = XDATCAR.readline().rstrip('\r\n').rstrip('\n')
        scaling_factor = float(XDATCAR.readline())
        lattice = np.zeros((3,3))
        for i in range(3):
            lattice[i] = np.array([float(j) for j in XDATCAR.readline().split()])
        list_ele = [j for j in XDATCAR.readline().split()]
        num_ele  = [int(j) for j in XDATCAR.readline().split()]
        tot_ele  = sum(num_ele)
             
        first_frame = []
        
        if not ('D' in XDATCAR.readline()):
            print('Not Direct coordinate!')
            raise ValueError()
            
        for _ in range(tot_ele):       
            first_frame.append([float(j) for j in XDATCAR.readline().split()])
        first_frame = np.array(first_frame)
        
        tot_nsw, last_frame = check_xdatcar_last_frame(os.path.join(path,"XDATCAR"),tot_ele)
        
        is_liquid = atoms_moved(first_frame,last_frame,lattice)
    
        ### create jobs
        ### if we pass the XDATCAR, it is hard to read code, because we need count
        XDATCAR.close()
        sel_nsw = create_job(path,tot_nsw,tot_ele,is_liquid,liquid_frames = tot_nsw//args.liquid_gap, solid_frames = tot_nsw//args.solid_gap)
        
#        dsq_jobs(path)
#        run_dsq(cwd,path,sel_nsw)
    

