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
from shared_functions import load_paths, check_outcar_done, check_incar

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


def create_job(path,tot_ele,sel_nsw):
    xdt2pos(path,sel_nsw,tot_ele,copybs = False)

def dsq_jobs(path,sel_nsw):
    """
    path: str
    sel_nsw: list
    """
    out = []
    target_folder = os.path.join(path,'recal')
    if len(sel_nsw)<1000:
        for i in sel_nsw:
            jobpath = os.path.join(target_folder,str(i+1))
            out.append('cd '+jobpath+'; ' + 'mpirun vasp_std' +'\n')
        fw   = open(os.path.join(target_folder,'out'),'w')
        fw.writelines(out)
        fw.close()        
    else:
        pass # consider if there are over 1000 frames

from subprocess import call
def run_dsq(cwd,path):
    target_folder = os.path.join(path,'recal')
    if 'out' in os.listdir(target_folder):
        os.chdir(target_folder)
        call("dsq --job-file out --mem-per-cpu 5g -t 08:00:00 -p scavenge --mail-type None --batch-file sub_script.sh\n",shell=True)
        call("sbatch sub_script.sh", shell=True)
        os.chdir(cwd)


def get_XDATCAR_info(path):
    
        XDATCAR = open(os.path.join(path,"XDATCAR"),'r')
        title   = XDATCAR.readline().rstrip('\r\n').rstrip('\n')
        scaling_factor = float(XDATCAR.readline())
        lattice = np.zeros((3,3))
        for i in range(3):
            lattice[i] = np.array([float(j) for j in XDATCAR.readline().split()])
        list_ele = [j for j in XDATCAR.readline().split()]
        num_ele  = [int(j) for j in XDATCAR.readline().split()]
        tot_ele  = sum(num_ele)
        
        if not ('D' in XDATCAR.readline()):
            print('Not Direct coordinate!')
            raise ValueError()

        XDATCAR.close()
        return scaling_factor, title, lattice, list_ele, num_ele, tot_ele


def read_outcar_sel_nsw(outcar):
    with open(outcar) as fp:
        for line in fp:
            if 'nsw_sel' in line:
                sel_nsw_outcar =  line.split('=')[1].replace('\n','').split()
                fp.close()
                break
    return [int(i) for i in sel_nsw_outcar]
    
def get_sel_nsw(relax_step):
    
    """
    output sel_nsw MUST be sorted!
    """
    if relax_step > 0:    
        sel1 = np.random.choice(range(20), 5, replace=False)
        sel2 = np.random.choice(range(20,100), 10, replace=False)
        sel3 = np.random.choice(range(100,200), 10, replace=False)
        sel4 = np.random.choice(range(200,500), 10, replace=False)
        sel5 = np.random.choice(range(500,1000), 10, replace=False)
        sel6 = np.random.choice(range(1000,2000), 10, replace=False)
        sel7 = np.random.choice(range(2000,5000), 10, replace=False)  #relax should not exceed 5000
        sel = np.concatenate((sel1, sel2,sel3,sel4,sel5,sel6,sel7))
        sel.sort()
        sel_nsw = sel[sel<relax_step]
    else:
        sel_nsw = []
    return sel_nsw
    
def build_recal(path, tot_ele, sel_nsw):
    
    if sel_nsw == []:
        print("No job, pass!")
    else:
        create_job(path,tot_ele,sel_nsw)
        dsq_jobs(path,sel_nsw)
        run_dsq(cwd,path)   
    
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--inputfile","-if",help="input files")

args   = parser.parse_args()

cwd    = os.getcwd()
if args.inputpath:
    print("Check files in {0}  ".format(args.inputpath))
    inputpath = args.inputpath
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd,'folders_pv')

if args.inputfile:
    print("Check files in {0}  ".format(args.inputfile))
    inputfile = args.inputfile
else:
    print("No folders point are provided. Use default value folders")
    inputfile = os.path.join(cwd,'inputs')
    
tmp = load_paths(inputpath)
paths = tmp[0]; relax_steps = tmp[1]; relax_paths = tmp[2]

call("dsq -h | head",shell=True)

for i in range(len(paths)):
    path = paths[i]
    print('###',path)
    relax_step = relax_steps[i]
    relax_path = relax_paths[i]
    scaling_factor, title, lattice, list_ele, num_ele, tot_ele = get_XDATCAR_info(path)
    if path == relax_path and os.path.exists(os.path.join(path,'recal')):
        sel_nsw_tmp = get_sel_nsw(relax_step)  # if relax_step==0 => [],
        outcar_sel_nsw =read_outcar_sel_nsw(os.path.join(os.path.join(path,'recal'),'OUTCAR'))
        sel_nsw = []
        for nsw in sel_nsw_tmp:
            if not (nsw in outcar_sel_nsw):
                sel_nsw.append(nsw)               
        build_recal(relax_path,tot_ele,sel_nsw)
    else:
        if not os.path.exists(os.path.join(relax_path,'recal')):
            os.mkdir(os.path.join(relax_path,'recal'))      
        sel_nsw = get_sel_nsw(relax_step)
        build_recal(relax_path,tot_ele,sel_nsw)
