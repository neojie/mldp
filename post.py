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
        print(target_path)
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
            print(coord)
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
            copy('inputs/'+target_incar,target_path)
            os.rename(os.path.join(target_path,target_incar),os.path.join(target_path,'INCAR'))
            copy('inputs/KPOINTS',target_path)
            copy('inputs/POTCAR',target_path)
#            if copybs:
#                copy('bs.sh',target_path)
#                copy('bs_job.sh',target_path)
    XDATCAR.close()

############################ read from file for paths  ############################
### get target path
def check_outcar_done(outcar):
    """
    check if outcar is done
    """
    txt  = reverse_readline(outcar)
    for i in range(int(1e4)):
        line = next(txt)
        if 'Voluntary context switches' in line:
#            print('OUTCAR is done')
            return True
    return False
kb = 8.617333262e-5
def check_incar(incar):
    """
    check if outcar is done
    """
    with open(incar) as incar_file:
        for line in incar_file:
            if 'SIGMA' in line:
               sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            if 'TEBEG' in line:
               tebeg = float(re.sub('[^0-9]', '', line.split('=')[1]))
    if abs(kb*tebeg - sigma)> 0.01:
        print(path,'SIGMA wrong')
        raise ValueError()
    else:
        sel_incar='INCAR_'+str(int(round(tebeg/1000)))+'k'
        return sel_incar
      
def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order
    read big file in reverse order is not trivial, library avaialable but not a optimum choice
    source: https://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python
    """
    
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment
            
#def read_path():
#    """
#    """
#    paths
#    return paths

def get_info(path,tot_ele):
    last_frame = []
    XDATCAR_reverse = reverse_readline(os.path.join(path,'XDATCAR'))
    for _ in range(tot_ele):
        line =  next(XDATCAR_reverse)
        last_frame.append([float(j) for j in line.split()])
    tot_nsw  = int(re.sub("[^0-9]", "",next(XDATCAR_reverse)))
    last_frame.reverse()
    return tot_nsw, np.array(last_frame)



def create_job(path,tot_nsw,tot_ele,is_liquid,liquid_frames = 500, solid_frames = 100,copybs = True):
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

def load_paths(inputpath):
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    fp.close()
    for i in range(len(folders_org)):
        if '#' in folders_org[i] or folders_org[i] == '\n':
            pass
        
        else:
            tmp = folders_org[i].replace('\n','')
            paths.append(tmp.replace('"','').replace(',','').replace('deepmd',''))
    return paths

def dsq_jobs(path,sel_nsw):
    out = []
    target_folder = os.path.join(path,'recal')
    if len(sel_nsw)<1000:
        jobpath = os.path.join(target_folder,str(sel_nsw+1))
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
    
def sub_job(path):  
    """
    """
    XDATCAR=open(os.path.join(path,"XDATCAR"),'r')
    
    
#paths = read_path()
#paths = ['/Users/jiedeng/Documents/tmp/jd848/Fe98Si10O20/melt4/r4-4000-FeSiC/']
import argparse
parser = argparse.ArgumentParser()
args   = parser.parse_args()
cwd    = os.getcwd()
parser.add_argument("--inputpath","-ip",help="input path file")
if args.inputpath:
    print("Check files in {0}  ".format(args.path))
    inputpath = args.inputpath
else:
    print("No folders point are provided. Use default value folders")
    inputpath = os.path.join(cwd,'folders')
    
#paths = ['/gpfs/loomis/project/kklee/jd848/pv+hf/4k/100g/r3']
### 


paths = load_paths(inputpath)
for path in paths:
    # get basic info

    if os.path.exists(os.path.join(path,'recal')):
#        path = os.path.join(path,'recal')
#        pass 
        print(path,'has recal folder')
    else:
        os.mkdir(os.path.join(path,'recal'))        
        XDATCAR = open(os.path.join(path,"XDATCAR"),'r')
        title   = XDATCAR.readline().rstrip('\r\n').rstrip('\n')
        scaling_factor = float(XDATCAR.readline())
        lattice = np.zeros((3,3))
        for i in range(3):
            lattice[i] = np.array([float(j) for j in XDATCAR.readline().split()])
        list_ele   = [j for j in XDATCAR.readline().split()]
        num_ele = [int(j) for j in XDATCAR.readline().split()]
        tot_ele = sum(num_ele)
             
        first_frame = []
        
        if not ('D' in XDATCAR.readline()):
            print('Not Direct coordinate!')
            raise ValueError()
            
        for _ in range(tot_ele):       
            first_frame.append([float(j) for j in XDATCAR.readline().split()])
        first_frame = np.array(first_frame)
        
        tot_nsw, last_frame = get_info(path,tot_ele)
        
        is_liquid = atoms_moved(first_frame,last_frame,lattice)
    
        ### create jobs
        ### if we pass the XDATCAR, it is hard to read code, because we need count
        XDATCAR.close()
        sel_nsw = create_job(path,tot_nsw,tot_ele,is_liquid,liquid_frames = tot_nsw//20, solid_frames = tot_nsw//40)
        dsq_jobs(path,sel_nsw)
        run_dsq(cwd,path)
    


