#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 23:36:41 2022

merge the soap similarity analyzed file of each element into one .xyz file

@author: jiedeng
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--output", "-o", type=str,default='merge.xyz',help="output file name")
parser.add_argument("--index", "-i",type=int,default = 10000,help="up to index i will the results be merged")
parser.add_argument('--ele',"-e", type=str,nargs="+",default = ['Mg','Si','O','H'],help="elements")
parser.add_argument('--file',"-f", type=str,nargs="+",
                    default = ['refmg-local.xyz','refsi-local.xyz',
                               'refo-local.xyz','refh-local.xyz'],help="path to file")
args   = parser.parse_args()

def modify_cell(cell):
    cs=cell.split()
    cs.insert(1,'0.0000');cs.insert(1,'0.0000');cs.insert(1,'0.0000')
    cs.insert(5,'0.0000');cs.insert(5,'0.0000');cs.insert(5,'0.0000')
    cs2 = '{0}=\"{1}\"'.format('Lattice'," ".join(map(str,cs)))
    return cs2

def comment(cs2):
    return cs2+" Properties=species:S:1:pos:R:3:k:R:1\n" # k is the similarity kernel

def merge(objects,num,ele,outname='out.xyz'):
    out = open(outname,'w')
    tot_ele_num = len(objects)
    atom_num = []
    for jj in range(num):
        try:
            # tot_atom_num = 0
            for ii in range(tot_ele_num):
                atom_num.append(int(objects[ii].readline().rstrip('\r\n').rstrip('\n')))
                cell = objects[ii].readline()
            tot_atom_num = sum(atom_num)
            out.writelines(str(tot_atom_num)+'\n')
            ct = comment(modify_cell(cell))
            out.writelines(ct)
            
            for kk in range(tot_ele_num):
                for ll in range(atom_num[kk]):
                    out.writelines(objects[kk].readline().replace('X',ele[kk]))
            atom_num = []
        except:
            print('read {0} < {1}'.format(jj,num))
            break
    for obj in objects:
        obj.close()
    out.close()

"""
f1 = open('/Users/jiedeng/GD/papers/pv7_h/partition/ppv-deng/ppv-NN-JD/refmg-local.xyz')
f2 = open('/Users/jiedeng/GD/papers/pv7_h/partition/ppv-deng/ppv-NN-JD/refsi-local.xyz')
f3 = open('/Users/jiedeng/GD/papers/pv7_h/partition/ppv-deng/ppv-NN-JD/refo-local.xyz')
f4 = open('/Users/jiedeng/GD/papers/pv7_h/partition/ppv-deng/ppv-NN-JD/refh-local.xyz')
merge([f1,f2,f3,f4],10,['Mg','Si','O','H'],'test.xyz')
"""
objects = []
for file in args.file:
    objects.append( open(file))

merge(objects,args.index,args.ele,args.output)
    
    