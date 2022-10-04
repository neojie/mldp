#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 07:28:21 2020

@author: jiedeng
"""

import os
#import re
#import shutil
import glob
import argparse
from subprocess import call

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str,help="input file")
parser.add_argument("--n_frame","-n",type=int,default= 1000, help="number of frames")
parser.add_argument("--tail","-t", default=True,action='store_false', help="tail of dump file? default YES")
parser.add_argument('--range',"-r", type=int,nargs="+",help="range of index starting from 0")

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
                
args   = parser.parse_args()

if args.file:
    file = args.file
else:
    print("No file supplied, search for *dump file")
    files  = glob.glob("*dump")
    file = files[0]
    print("Find {0}; analyze {1}".format(files,file))

# if args.tail:
txt  = reverse_readline(file)
# else:
#     txt=open(file)
    # txt  = fp.readline()
    
# lines = []
# count = 0

# for i in range(int(1e6)):   #somehow for some cases there is 'Voluntary context switches' in between?!
#     if args.reverse:
#         line = next(txt)
#         if 'TIMESTEP' in line:
#             break
#     else:
#         line = txt.readline()
#         if count ==2:
#             break
#         if 'TIMESTEP' in line:
#             count+=1            
    # lines.append(line)
for i in range(int(1e6)):   #somehow for some cases there is 'Voluntary context switches' in between?!
    line = next(txt)
    if 'TIMESTEP' in line:
        break
  


print('Total # of atoms is {0}'.format(i-8))

try:
    n = len(args.range)
    if n==1:
        print("from frame 0 to {}".format(args.range[0]))
        beg = 0 
        end = args.range[0]
    elif n==2:
        beg, end = min(args.range), max(args.range)
        print("from frame {} to {}".format(beg,end))
    else:
        print("ERROR")
        import sys
        sys.exit('length of range is WRONG')

    com = "sed -n '{0},{1}p;{2}q' tail_10.dump > range_{3}_{4}.dump".format((i+1)*beg+1,(i+1)*end,(i+1)*end+1,beg,end)
    # print(com)
    call(com,shell=True)
        
except:
    if args.tail:
        call("tail -n {0} {1} > tail_{2}.dump".format((i+1)*args.n_frame,file,args.n_frame), shell=True)
    else:
        call("head -n {0} {1} > head_{2}.dump".format((i+1)*args.n_frame,file,args.n_frame), shell=True)

