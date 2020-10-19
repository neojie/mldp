#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:17:57 2020
when Pressure is too high, VASP OUTPUT format for stress are incorrect!
@author: jiedeng
"""
import os

def modify(line):
    line_split = line.split()
    tmp = line_split[2]
    ll = len(tmp)
    p1 = tmp[:ll//3]
    p2 = tmp[ll//3:ll*2//3]
    p3 = tmp[ll*2//3:]
    ave = (float(p1) + float(p2) + float(p3))/3
    if abs(float(p1) - float(p2)) < ave and abs(float(p1) - float(p3)) < ave and abs(float(p2) - float(p3)) < ave:
        pass
    else:
        print(p1,p2,p3,"too different")
        raise ValueError()
    
    
#    line_split[2] = ' '.join([p1,p2,p3])
    
    return line.replace(tmp,' '.join([p1,p2,p3]))

outcar1 = 'OUTCAR'
outcar2 = 'OUTCAR2'
fp1 = open(outcar1)
fp2 = open(outcar2,'w')

for line in fp1:
#    print(line)
    if ('in kB' in line) and len(line.split())==6:
#        print(line)
        line2 = modify(line)
#        print(line2)
        fp2.write(line2)
        fp2.write('\n')
    else:
#        print(")
        fp2.write(line)
fp1.close()
fp2.close()
os.rename('OUTCAR','OUTCAR_org')
os.rename('OUTCAR2','OUTCAR')

#line='in kB  396394.48825395676.36574391590.20644  1062.09175   725.12530     7.41405'
#line_split = line.split()
# focus on 3rd element


