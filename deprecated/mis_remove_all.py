#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 12:43:21 2020
for i in *;do rm $i/POTCAR $i/XDATCAR $i/CONTCAR *out;done
@author: jiedeng
"""
import os
import argparse
import shutil
import glob
from shared_functions import merge, merge_sel, remove_recal_traj_files

parser = argparse.ArgumentParser()
parser.add_argument("--delete","-d",help="to delete")
parser.add_argument("--appenddelete","-ad",help="append to delete")
args   = parser.parse_args()

cwd    = os.getcwd()

if args.delete:
    print("To delete {0}  ".format(args.delete))
    delete = args.delete
    delete = delete.split()
else:
    print("No Files provided!")
    delete = ['PROCAR','CHGCAR','CHG', 'REPORT','WAVECAR','bsfoohpymelt', \
              'EIGENVAL','IBZKPT', 'PCDAT','pet.dat','sub_job.sh','DOSCAR','vasprun.xml',\
              'bsfoohpymelt-grace','dipole','inp.vis','XDATCAR.xzy','unwrap*']


# remove all db files
# mark train files
# remove deepmd file 
# remove traj/recal file => merge to OUTCAR
if args.appenddelete:
    for item in args.appenddelete.split():
        if not item in delete:
            delete.append(item)
# delete current folder, all sub folder, sub folder of the subfolder
def remove_file(path,wild_card):
    target =  os.path.join(path,wild_card)
    for file_to_remove in glob.glob(target):
        os.remove(file_to_remove)
        
def delete_file(path,to_delete):
    """
    """
    files = os.listdir(path)
    for file in to_delete:
        if file in files:
            os.remove(os.path.join(path,file))
    for file in files:
        file = os.path.join(path,file)
        if '.db' in file: # remove .db file
            os.remove(file)
        if 'train' in file:
            ask_for_delete.append(file)
        if os.path.isdir(file) and 'traj' in os.path.basename(file) and not os.path.exists(os.path.join(file,'OUTCAR')):
            print("### remove traj",file)
            merge_sel(file);remove_recal_traj_files(file)            
        if os.path.isdir(file) and 'recal' in os.path.basename(file) and not os.path.exists(os.path.join(file,'OUTCAR')):  ## for the recal folder only
            print("*** remove recal",file)
            merge_sel(file);remove_recal_traj_files(file);remove_file(file,'outcar')         
        if 'DG' in file:
            try:
                shutil.rmtree(file)
            except:
                print("-->",file)
            
directories  = [x[0] for x in os.walk(cwd)]
ask_for_delete = []
for directory in directories:
#    print(directory)
    if "deepmd" in directory and not('recal' in directory):
        if os.path.exists(directory):
            shutil.rmtree(directory)
    if os.path.exists(directory):
        delete_file(directory,delete)
print("###"*30)
print(ask_for_delete)
#delete_or_not  =  input("Delete these file (default is Y): ")
delete_or_not  = 'y'
if not (delete_or_not == 'N'):
    for path in ask_for_delete:
        if os.path.exists(path):
            try:
                shutil.rmtree(path) 
            except:
                os.remove(path)


