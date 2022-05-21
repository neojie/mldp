#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:26:56 2022

@author: jiedeng
"""

def cat(outfilename, *infilenames):
    with open(outfilename, 'w') as outfile:
        for infilename in infilenames:
            with open(infilename) as infile:
                for line in infile:
                    if line.strip():
                        outfile.write(line)
# cat('stat_6k1.py', 'stat_6k.py', 'stat_6k.py')

def cat2(outfilename, *infilenames):
    with open(outfilename, 'w') as outfile:
        for infilename in infilenames:
            with open(infilename) as infile:
                for line in infile:
                    if line.strip():
                        outfile.write(line)
# name='stat_{0}_{1}.txt'.format(int(ch[:,0][0]),int(ch[:,0][-1]))
# save_ch(ch,0,2.5,os.path.abspath('.'),name)

log = open('log.sub','r')


with open('sum_counts.txt', 'w') as outfile:
    for line in log:
        f1 = line.strip().replace('sh','txt').replace('stat','sum_counts')
        with open(f1) as infile:
            for line in infile:
                if line.strip():
                    outfile.write(line)

log = open('log.sub','r')
with open('sum_dimensions.txt', 'w') as outfile:
    for line in log:
        f1 = line.strip().replace('sh','txt').replace('stat','sum_dimensions')
        with open(f1) as infile:
            for line in infile:
                if line.strip():
                    outfile.write(line)