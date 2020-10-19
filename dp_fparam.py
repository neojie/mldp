#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 09:11:01 2020

make deep MD raw data
make temperature data
bs_fparam.py
bs_
read from file

@author: jiedeng
"""

import numpy as np
import glob
import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--version","-v",help="show program version",action= "store_true")
parser.add_argument("--path","-p",help="path file")

args = parser.parse_args()

if args.version:
    print("This is  version 0.1. xx")
    
cwd  = os.getcwd()

if args.version:
    print("This is  version 0.1. make deep MD raw data and temperature files")
    
if args.path:
    print("Check files in {0}  ".format(args.path))
    path = args.path
else:
    print("No folders point are provided. Use default value folders")
    path = os.path.join(cwd,'folders')

big_marker = '#'*100
small_marker = '-'*100

fp          = open(path,'r')
folders_org = fp.readlines()
folders     = []
fp.close()



for i in range(len(folders_org)):
    if '#' in folders_org[i] or folders_org[i] == '\n':
        pass
    else:
        tmp = folders_org[i].replace('\n','')
        folders.append(tmp.replace('"','').replace(',',''))

print('-'*40, "{0} folders selected are ".format(len(folders)),'-'*40)
for path in folders:
    print(path)
print(small_marker)




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
            
## memory saving way
def check_outcar(path):
    exist_outcar = True
    nsw = 0
    outcar_done = False
    if 'deepmd' in path: # if deepmd in path, remove deepmd
        path = path.split('deepmd')[0]
    outcar = os.path.join(path,'OUTCAR')
    if not os.path.exists(outcar):
        exist_outcar = False
    else:   
        txt  = reverse_readline(outcar)
        for i in range(int(1e4)):
            line = next(txt)
            if 'Voluntary context switches' in line:
                outcar_done = True
            if 'Iteration' in line:
                nsw = int(line.split('Iteration')[1].split('(')[0])
                break
    return exist_outcar, nsw, outcar_done

print('-'*40, 'Check the OUTCAR','-'*40)


nsws = []
all_outcar_good = True
for path in folders:
    exist_outcar, nsw, outcar_done = check_outcar(path)
    nsws.append(nsw)
    if not exist_outcar:
        print(path, 'no outcar')
#        raise FileNotFoundError()
        all_outcar_good = False
    if nsw < 100:
        print('{0} has {1} iteration only'.format(path,nsw))
#        raise ValueError()

#if not all_outcar_good:
#    raise FileNotFoundError()
print('-'*40, "All folders have OUTCAR with iter >= {0} ".format(100),'-'*40)

print('-'*40, 'check and build deepmd file','-'*40)


# check if deepmd folder exists, -> does -> check how many sets, each set size -> liquid only 2000 at most for train,1000 for test
# if more than 3 sets, change others into xsets

# no deepmd, form deepmd
from dpdata import LabeledSystem
def build_deepmd(path,nsw):
    ls=LabeledSystem(os.path.join(path, 'OUTCAR'),fmt='outcar')
    deepmd = os.path.join(path,'deepmd')
    if nsw <= 2000: # we know nsw must > 100
        set_size = nsw//2
    if nsw > 2000:
        set_size = 1000
    ls.to_deepmd_npy(deepmd,set_size=set_size)
    if nsw>3000:
        check_sets(deepmd)

def check_sets(deepmd):
    sets=glob.glob(deepmd+ "/set*")
    if len(sets) >3 and np.load(sets[0]+'/energy.npy').size>999:
        print(deepmd)
        print("{0} sets exist".format(len(sets)))
    for i in range(3,len(sets)):
        seti  = os.path.join(deepmd,'set.00{0}'.format(i))
        xseti = os.path.join(deepmd,'xset.00{0}'.format(i))
        # move set.002+ to xset.002
        os.rename(seti,os.path.join(deepmd,xseti))

kb = 8.617333262e-5
def build_fparam(path): # deepmd/..
    incar = os.path.join(path,'INCAR')
    
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
        deepmd = os.path.join(path,'deepmd')
        sets=glob.glob(deepmd+"/set*")
        for seti in sets:
            energy=np.load(seti+'/energy.npy')
            size = energy.size
            all_te = np.ones(size)*sigma
            np.save( os.path.join(seti,'fparam.npy'), all_te)
  
def check_deepmd(path,nsw):

    if 'deepmd' in path or os.path.exists(os.path.join(path,'deepmd')): # if deepmd in path, remove deepmd
        print("deepmd foler already exist")
        check_sets(path)
        if 'fparam.npy' in os.listdir(os.path.join(path,'set.000')):
            pass
        else:
            build_fparam(path.split('deepmd')[0])
    else:
        build_deepmd(path,nsw)
        build_fparam(path)
        
for i in range(len(folders)):
    path = folders[i]
    print(path)
    check_deepmd(path,nsws[i])
    
print('deepmds are done')
print(big_marker)


#for i in range(len(folders)):
    
    

