#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 21:58:14 2022

@author: jiedeng
"""
import argparse
from subprocess import call
import os
parser = argparse.ArgumentParser(description="build TI files\
                                 la1 must exist")
parser.add_argument("--nsw","-n",type=int,default=5000, help="nsw")
args = parser.parse_args()

scalees = [1.0 ,0.7179229, 0.3192687, 0.0808200, 0.009658, 0.0003546, 0.0000011]

pwd = os.path.abspath(os.curdir)
com = "mkdir la{2..7};for i in la{2..7};do cp la1/INCAR la1/CONTCAR la1/KPOINTS la1/POTCAR $i;mv $i/CONTCAR $i/POSCAR;done"
call(com,shell=True)

def key_val_sub(key,val):
    com =  '/{0}/s/\([0-9].\+\)/{1}/g'.format(key,val)
    com = "sed -i -e '{0} ' {1}".format(com,'INCAR')
    print(com)
    call(com,shell=True) # sed to replace numbers after SIGMA, 
    
for i in range(1,7):
    print(i)
    dire = 'la{0}'.format(i+1)
    os.chdir(dire)
    key_val_sub('SCALEE',scalees[i])
    key_val_sub('NSW',args.nsw)
    os.chdir(pwd)

print("**to submit jobs**")
print("for i in la{2..7};do cd $i;qsubvaspti;cd ..;done")

