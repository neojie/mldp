#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:36:59 2022

@author: jiedeng
"""


## generate index file, absolute index from 0
import numpy as np
import os, re
from shutil import copyfile
import fileinput
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--begin","-b",type = int, default = 500, help="begin index, default: 500")
parser.add_argument("--end","-e",type = int, default = 4000, help="end index, default: 4000")
parser.add_argument("--step","-s",type = int, default = 30, help="step, default: 30")
args   = parser.parse_args()


index = np.array(range(args.begin,args.end,args.step))
if os.path.exists('index'):
    print('***step 1: index file exists, skip')
else:
    print('---step 1: generate index file')
    np.savetxt('index',index,fmt='%d')
## use extract deepmd to build deepmd based on index file
from subprocess import call
if os.path.exists('deepmd'):
    print('***step 2: deepmd file exists, skip')
else:
    print('---step 2: generate deepmd file') 
    call('python ~/script/mldp/extract_deepmd.py -ttr 10000000 -id index -st',shell=True)

## build inputs file

def key_word_replace(fname='INCAR',key_word_search='ENCUT',replacement_text='600'):
    """
    replace all the number for key word
    key_word_search : key word to search
    the problem come from 400.00 and 400 is different
    """
    tx_org = 'no original word found'
    tx_rep = 'no replacement done'    
    ## wild card is used here
#    for f_ele in os.listdir('.'):
#        if fnmatch.fnmatch(f_ele, fname):
#            f_full = f_ele
    f_full = os.path.abspath(fname)
    if os.path.isfile(f_full):
        with fileinput.FileInput(f_full, inplace=True, backup='-bak') as file:
            for line in file:
                if key_word_search in line:
                    
                    tx_org = line
                    #line=line.split()
                    #line_split = =line.split()
                    num = re.findall('\d*\.?\d+',line)
                    text_to_replace = ' '.join(map(str, num))
                    print(line.replace(text_to_replace, replacement_text))
                    tx_rep = replacement_text
                else:
                    print(line, end='')
#            os.rename(f_full+'-bak','old'+f_full)
    else:
        print('No'+ f_full+ 'found')
    return tx_org, tx_rep

print("Build iputs")

if os.path.exists('inputs'):
    print('***step 3: inputs for recal exist, skip')
else:
    print('---step 3: generate inputs for recal')  
    os.mkdir('inputs')
    copyfile('INCAR', 'inputs/INCAR')
    copyfile('KPOINTS', 'inputs/KPOINTS')
    copyfile('POTCAR', 'inputs/POTCAR')
    tx_org, tx_rep = key_word_replace(fname='inputs/INCAR',key_word_search='NSW', replacement_text=str(0))
    tx_org, tx_rep = key_word_replace(fname='inputs/INCAR',key_word_search='ENCUT', replacement_text=str(800))
    call("sed -i -e 's/1/{0}/g' inputs/KPOINTS".format(2), shell=True)

## recal using deepmd_recal
if os.path.exists('recal/1') or os.path.exists('recal/OUTCAR') or os.path.exists('recal/toten.dat'):
    print('step 4: recaled results exist, skip')
else:
    call('python ~/script/mldp/recal_dpdata.py -d deepmd',shell=True)
    import sys
    sys.exit()
# post recal analis
print('---step 5: check if all runs are done, inside the recal folder')
print('python ~/script/mldp/post_recal.py')
print('---step 6: check the convergence ')
print('python ~/script/mldp/check_nbands_nelm.py -ip all -v')
print('---step 7: merge all the runs to a single OUTCAR')
print('python ~/script/mldp/merge_out.py -o OUTCAR -r y')
print('---step 8: extract toten')
print('pet -l toten')
print('extract the results')

## analyze results, toten.dat in and recal/totend.dat
if os.path.exists('recal/OUTCAR') and os.path.exists('recal/toten.dat') and os.path.exists('toten.dat'):
    import matplotlib.pyplot as plt

    print('---step 9: thermal dynamic pertubation, df = target - ref')
    e_org   = np.loadtxt('./toten.dat')
    e_recal = np.loadtxt('./recal/toten.dat')
    idx_file = 'index'
    idx = np.loadtxt(idx_file, dtype=int)
    
    #
    t = float(input('input tempererature: '))
    boltz = 8.617e-5
    beta = 1/t/boltz
    du  = e_recal - e_org[idx]
    df1=-1/beta*np.log(np.average(np.exp(-du*beta)))
    df2=np.average(du) - beta/2*np.average(du**2) + beta/2*np.average(du)**2
    print("du_mean, ev/system",np.average(du),"df1 ev/system",df1,"df2 ev/system",df2)
    
    df1_step, df2_step = np.zeros(idx.shape), np.zeros(idx.shape)
    for i in range(len(idx)):
        du = e_recal[:i+1] - e_org[idx][:i+1]
        df2_step[i] = du.mean() - beta/2*np.mean(du**2) + beta/2*(du.mean()**2)
        df1_step[i] = -1/beta*np.log(np.average(np.exp(-du*beta)))

    plt.figure()
    fig,ax = plt.subplots(1,2,figsize=(11,4),sharey=False)
    ax[0].plot(e_org[idx],label='org toten')
    ax[0].plot(e_recal,label='recal toten')
    ax[0].legend()
    ax[1].plot(range(len(idx)),df1_step,label='df1, exp way,eqn 190 Dorner')
    ax[1].plot(range(len(idx)),df2_step,label='df2, exp way,eqn 191 Dorner')
    ax[1].set_xlabel('number of snapshots');ax[1].set_ylabel('dF')
    fig.show()
    fig.savefig('tpt.pdf',bbox_inches='tight')


else:
    print('***step 9: thermodynamic perturbation, skip due to lack of files')