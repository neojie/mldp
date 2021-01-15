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
parser.add_argument("--timestep", "-t",type=float, default=1,  help=" timestep in fs, default 1fs")
parser.add_argument("--scale", "-s",type=float, default=1,  help=" scale to SI unit, check the log file for this value, default 1")

args = parser.parse_args()

dat=read_file(args.input_file,args.num)

dt = dat[:,1]*args.timestep/1e3 # in ps
JxJx = dat[:,3] # autocorr of heat flux in x direction
JyJy = dat[:,4]
JzJz = dat[:,5]
JJ = (JxJx + JyJy + JzJz)/3
JJ_JJ0 = JJ/JJ[0]

cumsum_JxJx = np.cumsum(JxJx)*args.scale
cumsum_JyJy = np.cumsum(JyJy)*args.scale
cumsum_JzJz = np.cumsum(JzJz)*args.scale

cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3

print('integrations of x,y,z are : ', cumsum_JxJx[-1], cumsum_JyJy[-1], cumsum_JzJz[-1])
print('This number should be consistent with log file value: ', cumsum_JxJx[-1], cumsum_JyJy[-1], cumsum_JzJz[-1])
print(' kappa is {0} (W m-1 K-1): ', cumsum_JJ[-1])

fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
ax[0].plot(dt,JxJx,label='x') # 2ps is enough, interesting
ax[0].plot(dt,JyJy,label='y')
ax[0].plot(dt,JzJz,label='z')
ax[0].plot(dt,JJ,label='average')
ax[0].set_xlabel('dt (ps)')
ax[0].ylabel('autocorrelation')
ax[0].grid(True)
ax[0].legend()

ax[1].plot(dt,cumsum_JxJx,label='x') # 
ax[1].plot(dt,cumsum_JyJy,label='y')
ax[1].plot(dt,cumsum_JzJz,label='z')
ax[1].plot(dt, cumsum_JJ,label='average')

ax[1].set_xlabel('dt (ps)')
ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
ax[1].grid(True)
ax[1].legend()
plt.show()
#plt.figure()
#plt.plot(dt,JJ_JJ0)
#plt.show()