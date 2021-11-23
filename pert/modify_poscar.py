#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 12:48:28 2021

p0-x: append all newatom position in the end also change the title line

all: append all newatom position in the end also change the title line

After making the files

for i in {0..9};do echo $i;mkdir $i;cp p$i.vasp $i/POSCAR;cp /u/project/ESS/lstixrud/jd848/metad/pvh/inputs/inputs_4000/* $i;cd $i;qsub sub_vasp.sh;cd ..;done
for i in {0..8};do echo $i;mkdir $i;cp p$i.vasp $i/POSCAR;cp /u/project/ESS/lstixrud/jd848/metad/pvh/inputs/inputs_4000/* $i;cd $i;qsub sub_vasp.sh;cd ..;done

@author: jiedeng
"""
print('--'*40)
print('prepare pristine cell as p0.vasp')
print('assign two indexs of atoms, starting from 0')
print('specify number of poscar file to be generated')
print('--'*40)
import numpy as np
import copy
import argparse

parser = argparse.ArgumentParser(description="modify poscar")
parser.add_argument("-i", type=int, nargs="+", help="index, the 1st approaches the 2nd")
parser.add_argument("-num", type=int, help="num")

args = parser.parse_args()

poscar='p0.vasp'
datain = open(poscar, 'r')
idx1 = args.i[0] # python index
idx2 = args.i[1] #
num = args.num # to get num atoms in between atom 1 and 2

dat = datain.readlines()
newdat = copy.copy(dat)
eles = dat[5].split()
nums = dat[6].split()
nums_tmp = np.array(nums).astype(int)
nums_acc = np.cumsum(nums_tmp)

def judge_ele(idx):
    for i in range(len(nums_acc)):
        if idx+1 <= nums_acc[i]:
            return eles[i]

ele1=judge_ele(idx1)
ele2=judge_ele(idx2)
            
lastlineidx = int(np.array(nums).astype(int).sum())+8

atom1 = dat[8+idx1]
atom2 = dat[8+idx2]

atom1 = np.array(atom1.split()).astype(float)
atom2 = np.array(atom2.split()).astype(float)

vec = atom2 - atom1

dv  = vec/(num +1)  # 

newatoms = []

def savedat(newdat,label):
    fp = open(label,'w')    
    for line in newdat:
        fp.write(line)
    fp.close()
    
for i in range(num):# two atom already known
    idx = i+1
    newatom = atom1 + dv*idx
    newatom_str = [str(jj) for jj in newatom]
    newatom_str = '    '.join(newatom_str)
    newatom_str = '    ' + newatom_str + '\n'
    newatoms.append(newatom_str)
    newdat[8+idx1] = newatom_str  
    label = 'p'+str(idx)+'.vasp'
    savedat(newdat, label)
    
for j in range(len(newatoms)):
    idx = lastlineidx + j
    dat.insert(idx,newatoms[j])
    
neweles = copy.copy(eles)
neweles.append('Xx')
neweles = '    '.join(neweles)
neweles = '    ' + neweles + '\n'

newnums = copy.copy(nums)
newnums.append(str(len(newatoms)))
newnums = '    '.join(newnums)
newnums = '    ' + newnums + '\n'

dat[5] = neweles
dat[6] = newnums

savedat(dat,'all.vasp')

### write log

log = [
       '# 0.vasp is pristine cell; all.vasp is trajectory \n',
       '# np.loadtxt(readme).astype(int) extract the (python, starting from 0) index of atoms of interests \n',
       '# {0} approaches {1} \n'.format(ele1, ele2),
       "# inspect file by  for i in *.vasp;do echo $i;sed '{0},{1}!d' $i;done \n".format(9+idx1,9+idx1),
        str(idx1) + '\n',
        str(idx2) + '\n',
       ]
savedat(log,'log')