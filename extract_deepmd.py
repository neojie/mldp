#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:51:58 2020
1- support dump data

@author: jiedeng
"""

import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("--inputpath","-ip",help="input path file")
parser.add_argument("--train_test_ratio","-ttr",type = float,default=3, help="input path file")
#parser.add_argument("--OUTCAR","-o",type = str, default = 'OUTCAR', help="OUTCAR name")
parser.add_argument("--file","-f",type = str, default = 'OUTCAR', help="file name, default: OUTCAR")
parser.add_argument("--temp","-t",type = float, help="temperature in K, only needed if file does not contain temperature info")
parser.add_argument("--format","-fmt",type = str, default = 'outcar', help="format, e.g., outcar, dump supported by DPDATA,")
parser.add_argument("--deepmd","-d",type = str, default = 'deepmd', help="deepmd folder name")
parser.add_argument("--vaspidx","-vid",type = str, help="idx file, vasp idx, idx[0] >= 1")
parser.add_argument("--idx","-id",type = str, help="idx file, idx[0] >= 0")
parser.add_argument('--savetest',"-st", default=True, action='store_false',help="Default: save test as set.001? ")
parser.add_argument('--force_limit',"-fl", type=float,nargs="+",help="force limit max and min, order does not matter.")
parser.add_argument('--exclude',"-ex", type=int,nargs="+",help="manually exclude indexs")
parser.add_argument('--elements',"-e", type=str,nargs="+",help="impose element list")

args   = parser.parse_args()


import os
from dpdata import LabeledSystem
from dpdata import System
import glob
import numpy as np
import shutil


def extract_sigma_outcar(outcar):

    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    outcar_file = open(outcar)
    for _ in range(1000000):
        line = outcar_file.readline()
        if 'SIGMA' in line:
            sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            outcar_file.close()
            return sigma
        
def build_fparam(path, sigma, deepmd): # deepmd/..
#    
    sets=glob.glob(deepmd+"/set*")
    for seti in sets:
        nframe=np.load(seti+'/coord.npy').shape[0]
        print("set {0} : {1}".format(seti,nframe))
        all_te = np.ones(nframe)*sigma
        np.save( os.path.join(seti,'fparam.npy'), all_te)


def build_deepmd(path,outcar,deepmd):
    build_deepmd_frames(path,outcar,deepmd)
    
    # get temperature!
    static= False
    if args.temp is None:       
        if args.format=='outcar':
            sigma = extract_sigma_outcar(outcar)
        else:
            static= True
    else:
        kb = 8.617333262e-5
        sigma = args.temp*kb
    
    if not static:
        build_fparam(path,sigma, deepmd)
    else:
        print("No temperature info")

def update_element_files_with_defaults(type_map_dir='.', type_file_dir='.', new_elements_list=''):
    """
    Update the type_map.raw and type.raw files based on the new list of elements with default file paths.

    Parameters:
    type_map_dir (str): Directory path for type_map.raw. Default is current directory '.'.
    type_file_dir (str): Directory path for type.raw. Default is current directory '.'.
    new_elements_list (str): New list of elements (space-separated).
    """
    type_map_file_path = os.path.join(type_map_dir, 'type_map.raw')
    type_file_path = os.path.join(type_file_dir, 'type.raw')
    new_elements = new_elements_list#.split()

    # Read the current contents of the files
    with open(type_map_file_path, 'r') as file:
        current_type_map = file.read().splitlines()

    with open(type_file_path, 'r') as file:
        current_type = file.read().splitlines()

    # Store the original element list for summary
    original_element_list = current_type_map.copy()

    # Check if the new elements list is different from the current one
    if set(new_elements) != set(current_type_map):
        print("!!! Warning: Element list will be changed. !!!")

        # Update type_map.raw
        with open(type_map_file_path, 'w') as file:
            for element in new_elements:
                file.write(element + '\n')

        # Create a dictionary to map the new element indices
        element_index_map = {element: str(index) for index, element in enumerate(new_elements)}

        # Update type.raw with elements aligned in a single column
        updated_type = [element_index_map[current_type_map[int(atom)]] for atom in current_type]
        with open(type_file_path, 'w') as file:
            file.write('\n'.join(updated_type) + '\n')

        print(f"!!! Element list changed from {original_element_list} to {new_elements}. !!!")
    else:
        print("No changes needed. The element list is the same.")



def build_deepmd_frames(path,outcar,deepmd):
    """
    sub_ls = ls.sub_system(idx)
    
    """
    try:
        ls = LabeledSystem(outcar,fmt=args.format)
    except:
        ls = System(outcar,fmt=args.format)
        
    if args.exclude:
        oldsize = len(ls)
        idx_new = [i for i in range(len(ls)) if i not in args.exclude]
        ls = ls.sub_system(idx_new)
        newsize = len(ls)
        print('{0}/{1} is selected'.format(newsize,oldsize))
        
    if args.force_limit:
        fmin = min(args.force_limit)
        fmax = max(args.force_limit)
        print("force limit imposed, force in between {0}, {1}".format(fmin, fmax))
        idx_new = []
        exclude = []
        for i in range(len(ls)):
            forces = ls[i].data['forces']
            if forces.min() >= fmin and forces.max() <= fmax:
                idx_new.append(i)
            else:
                exclude.append(i)
        print('excluded frames', exclude)
        print('{0} / {1} is selected'.format(len(idx_new),len(ls)))
        ls = ls.sub_system(idx_new)

    if args.idx:
        print("index file provided")
        idx = np.loadtxt(args.idx).astype(int)
    elif (not args.idx) and  args.vaspidx:
        print("vasp index file provided")
        vaspidx=np.loadtxt(args.vaspidx)
        fp  = open(outcar)
        fp.readline()
        nsw_sel = fp.readline()
        if 'nsw_sel' in nsw_sel:
            print('file generated by merge_out.py')
            tmp     = nsw_sel.split('=')[1].strip().split(' ')
            nsw_sel = [int(tmp_idx) for tmp_idx in tmp]
            idx = []    
            for i in range(len(nsw_sel)):
                if nsw_sel[i] in vaspidx:
                    idx.append(i)
        else:
            print('OUTCAR file generated by VASP')
            idx = vaspidx - 1
    else:
        print("split train and test by ratio {0} : {1}".format(args.train_test_ratio,1))
        train_size = round(len(ls)*(args.train_test_ratio)/(args.train_test_ratio+1))
        idx = np.random.choice(range(len(ls)), train_size, replace=False)
        idx.sort()
        
    idx2 = [i for i in range(len(ls)) if i not in idx] # test
    ls2  = ls.sub_system(idx2) # test
    ls   = ls.sub_system(idx)
        
    deepmd = os.path.join(path,deepmd)
    
    ls.to_deepmd_npy(deepmd,set_size=1000000) # give a *large* value, default is 5000
    if len(ls2) == 0:
        print('test set has no data')
    elif args.savetest and len(ls2)>0:
        ls2.to_deepmd_npy('test_tmp',set_size=1000000)
        shutil.copytree('test_tmp/set.000',os.path.join(deepmd,'set.001'))
        shutil.rmtree('test_tmp')
    
    ## impose element list
    ## e.g.: OUTCAR only contains FeW, but the list in json is MgOFeW
    if args.elements:
        update_element_files_with_defaults( deepmd, deepmd, new_elements_list=args.elements)
    
if __name__ == "__main__":
    cwd    = os.getcwd()
    if os.path.exists(os.path.join(cwd,args.deepmd)):
        print("deepmd foler already exist=> skip")
    else:
        print("Build {0}".format(args.deepmd))
        build_deepmd(cwd,os.path.join(cwd,args.file), os.path.join(cwd,args.deepmd))
        if args.format =='dump':
            print("!!!** May want to modify the {0}/type_map.raw **!!!".format(args.deepmd))

