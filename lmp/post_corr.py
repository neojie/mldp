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
parser.add_argument("--timestep", "-ts",type=float, default=1,  help=" timestep in fs, default 1fs")
parser.add_argument("--scale", "-s",type=float,  help=" scale to SI unit, check the log file for this value, default 1")
parser.add_argument("--average", "-a",nargs="+",type=int, help=" step window average the thermal conductivity")
parser.add_argument("--temperature", "-t",type=float,help='temperature in K')
parser.add_argument("--volume", "-v",type=float,help='volume in A3')
parser.add_argument("--sample_rate", "-sr",type=float,help='sample rate')

args = parser.parse_args()

# convert from LAMMPS real units to SI
kB = 1.3806504e-23   # [J/K] Boltzmann
ev2j = 1.60218e-19
A2m = 1.0e-10
ps2s = 1.0e-12
convert = ev2j*ev2j/ps2s/A2m

if args.scale:
    scale = args.scale
    metal2SIps = scale/10/0.001 
    print('metal2SIps = scale/10/0.001 is assumed')
else:
    try:
        print("--"*40)
        print("  Unit conversion scale is not provided. Try parsing from inputs")
        if args.temperature:
            T = args.temperature
        else:
            import glob
            infile = glob.glob('in*')[0]
            print('  Find ', infile)
            print("  ?? temperature not provided, parse from in file") 
#            try:
            fp = open(infile)
            ins = fp.readlines()
            for line in ins:
                line=line.split('#')[0].split()
                if len(line)>0 and line[0] ==   'variable' and line[1] == 'T' and line[2] == 'equal':
                    T = float(line[3])
                    break
#            except:
#                print('No T found in ', infile)                   
            print('  ** T = ', T)
        if args.sample_rate:
            sr = args.sample_rate
        else:
            import glob
            infile = glob.glob('in*')[0]
            print('  Find ', infile)
            print(" ?? sample rate not provided, parse from in file") 
#            try:
            fp = open(infile)
            ins = fp.readlines()
            for line in ins:
                line=line.split('#')[0].split()
                if len(line)>0 and line[0] ==   'variable' and line[1] == 's' and line[2] == 'equal':
                    sr = float(line[3])
                    break
#            except:
#                print('No sample rate found in ', infile)   
            print(' ** sample rate = ', sr)
                
        if args.volume:
            V = args.volume
        else:
            infile = 'conf.lmp'
            fp = open(infile)
            ins = fp.readlines()
            print("  ?? vol not provided, parse from in conf.lmp") 
#            try:
            xhi = False
            yhi = False
            zhi = False
            for line in ins:
                if 'xlo' in line:
                    xlo = float(line.split()[0])
                    xhi = float(line.split()[1])
                if 'ylo' in line:
                    ylo = float(line.split()[0])
                    yhi = float(line.split()[1])
                if 'zlo' in line:
                    zlo = float(line.split()[0])
                    zhi = float(line.split()[1])
                if xhi and yhi and zhi:
                    break
                
            V = (xhi - xlo)*(yhi - ylo)*(zhi - zlo)
            print(' ** volume = ', V)

#            except:
#                print('volume parse error') 
            
        scale = convert/kB/T/T/V*sr*args.timestep/1e3
        metal2SIps = convert/kB/T/T/V
        print('  scale = ',scale)
        print("--"*40)
    except:
        raise ValueError('scale problem!')

dat=read_file(args.input_file,args.num)

dt = dat[:,1]*args.timestep/1e3 # in ps
JxJx = dat[:,3] # autocorr of heat flux in x direction
JyJy = dat[:,4]
JzJz = dat[:,5]
JJ = (JxJx + JyJy + JzJz)/3
JJ_JJ0 = JJ/JJ[0]

#cumsum_JxJx = np.cumsum(JxJx)*args.scale
#cumsum_JyJy = np.cumsum(JyJy)*args.scale
#cumsum_JzJz = np.cumsum(JzJz)*args.scale
#
#cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3
import scipy.integrate
cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)

cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3


print('integrations of x,y,z are : ', cumsum_JxJx[-1], cumsum_JyJy[-1], cumsum_JzJz[-1])
print('This number should be consistent with log file value: ')
print(' Last step kappa is {0} (W m-1 K-1): '.format( cumsum_JJ[-1]))

if args.average:
    window = args.average
    window_mean = cumsum_JJ[range(window[0],window[1])].mean()
    print('mean kappa within {0} ps - {1} ps is {2} (W m-1 K-1): '.format(
            dat[window[0],1]*args.timestep/1e3, 
            dat[window[1],1]*args.timestep/1e3, 
            window_mean))

 # s to ps

fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
ax[0].plot(dt,JxJx*metal2SIps,label='x') # 2ps is enough, interesting
ax[0].plot(dt,JyJy*metal2SIps,label='y')
ax[0].plot(dt,JzJz*metal2SIps,label='z')
ax[0].plot(dt,JJ,label='average')
ax[0].set_xlabel('dt (ps)')
ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
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

fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
ax[0].plot(dt,JxJx*metal2SIps,label='x') # 
ax[0].plot(dt,JyJy*metal2SIps,label='y')
ax[0].plot(dt,JzJz*metal2SIps,label='z')
ax[0].plot(dt,JJ,label='average')
ax[0].set_xlabel('dt (ps)')
ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
ax[0].grid(True)
ax[0].legend()
ax[0].set_xscale('log')

ax[1].plot(dt,cumsum_JxJx,label='x') # 
ax[1].plot(dt,cumsum_JyJy,label='y')
ax[1].plot(dt,cumsum_JzJz,label='z')
ax[1].plot(dt, cumsum_JJ,label='average')

ax[1].set_xlabel('dt (ps)')
ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
ax[1].grid(True)
ax[1].legend()
ax[1].set_xscale('log')
plt.show()

#plt.figure()
#plt.plot(dt,JJ_JJ0)
#plt.show()