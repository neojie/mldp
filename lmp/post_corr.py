#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 22:26:44 2021
analyze output of ave/correlate
@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
#from util import reverse_readline
import os 
import argparse

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
            
def read_file(file,count):
    out = []
    txt  = reverse_readline(file)
    for i in range(count):
        line = next(txt)
        line1 =[float(i) for i in line.split()]
        out.append(line1)
    return np.array(out[::-1])

parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
parser.add_argument("--input_file", "-i",type=str, default="J0Jt.dat",  help="ave/correlate output file")
parser.add_argument("--num", "-n",type=int, default=200,  help=" Nrepeat in ave/correlate Nevery Nrepeat Nfreq")
args = parser.parse_args()

print('timestep = 1s is assumed')
dat=read_file(args.input_file,args.num)

dt = dat[:,1]/1e3 # in ps
JxJx = dat[:,3] # autocorr of heat flux in x direction
JyJy = dat[:,4]
JzJz = dat[:,5]
plt.figure()
plt.plot(dt,JxJx,label='x') # 2ps is enough, interesting
plt.plot(dt,JyJy,label='y')
plt.plot(dt,JzJz,label='z')
JJ = (JxJx + JyJy + JzJz)/3
JJ_JJ0 = JJ/JJ[0]
plt.plot(dt,JJ,label='average')
plt.legend()
plt.grid(True)
plt.xlabel('dt')
plt.xlabel('correlation time (ps)')
plt.show()

cumsum_JxJx = np.cumsum(JxJx)
cumsum_JyJy = np.cumsum(JyJy)
cumsum_JzJz = np.cumsum(JzJz)
plt.figure()
plt.plot(dt,cumsum_JxJx,label='x') # 
plt.plot(dt,cumsum_JyJy,label='y')
plt.plot(dt,cumsum_JzJz,label='z')
print('integrations of x,y,z are : ', np.trapz(JxJx), np.trapz(JyJy), np.trapz(JzJz))
plt.grid(True)
#plt.figure()
#plt.plot(dt,JJ_JJ0)
#plt.show()