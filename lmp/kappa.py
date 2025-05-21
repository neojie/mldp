#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 11:02:36 2021
unit 
https://lammps.sandia.gov/doc/compute_heat_flux.html

directly analyze the log.properties file with Step, Jx, Jy, Jz 

On Wed Jan  6 22:26:44 2021
jiedeng: analyze output of ave/correlate

On Wed May 21, 2025
Yihang Peng: Improve the visualization

@author: jiedeng, Yihang Peng
"""

import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.integrate
import pickle as pkl
from matplotlib import rcParams

rcParams['font.size'] = '16'

parser = argparse.ArgumentParser(description="Calculate thermal conductivity from heat current autocorrelation function")
parser.add_argument("--input_file", "-i", type=str, default="jh.dat",  help="jh output file")
parser.add_argument("--timestep", "-ts", type=float, help="Timestep in fs, search in in.* if not provided")
parser.add_argument("--total_time", "-tt", type=float, default=100, help="total time for plotting in ps, default 100 ps")
parser.add_argument("--temperature", "-t", type=float,help='Temperature in K, search in Temp.txt and then in.* if not provided')
parser.add_argument("--volume", "-v", type=float,help='Volume in A^3, search in eq.final.lmp if not provided')
parser.add_argument("--configeration", "-conf", type=str, default='eq.final.lmp', help="read volume from this file")
parser.add_argument("--avg1", "-a1", type=float, default=0.1, help="average kappa between avg1 and avg2")
parser.add_argument("--avg2", "-a2", type=float, default=1, help="average kappa between avg1 and avg2")
parser.add_argument("--outfile", "-o", type=str,default='kappa', help="out file name")
parser.add_argument("--save_pkl", "-pkl", action='store_true', help="save data in pickle format")

args = parser.parse_args()

def autocorr(a):
    b=np.concatenate((a,np.zeros(len(a))),axis=0)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)//2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d

if args.temperature:
    T = args.temperature
elif os.path.exists('Temp.txt'):
    print("  ?? temperature not provided, parse from Temp.txt") 
    with open("Temp.txt", "r") as file:
        T = float(file.read())
        print('  ** T = ', T)
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? temperature not provided, parse from in file. Warning: The results may be not accurate!!!!") 
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

if args.timestep:
    timestep = args.timestep*1e-3#1e-3    # in ps this is different from post_corr
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? timestep not provided, parse from in file") 
    fp = open(infile)
    ins = fp.readlines()
    for line in ins:
        line=line.split('#')[0].split()
        if len(line)>0 and line[0] ==   'variable' and line[1] == 'dt' and line[2] == 'equal':
            timestep = float(line[3])
            break
#            except:
#                print('No T found in ', infile)                   
    print('  ** ts = ', timestep)

if args.volume:
    V = args.volume
else:
    infile = args.configeration
    fp = open(infile)
    ins = fp.readlines()
    print("  ?? vol not provided, parse from the configration file") 
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

# convert from LAMMPS real units to SI
kB   = 1.3806504e-23   # [J/K] Boltzmann
ev2j = 1.60218e-19
A2m  = 1.0e-10
ps2s = 1.0e-12
convert     = ev2j*ev2j/ps2s/A2m
scale = convert/kB/T/T/V*timestep
metal2SIps = convert/kB/T/T/V

file = args.input_file

J=np.loadtxt(file) # heat current J is saved, but heat flux (energy x velocity) autocorrelation is saved in J0Jt, heat flux/V = heat flux
J = J*V
Jx, Jy, Jz = J[:,1], J[:,2], J[:,3]
t = np.arange(len(Jx)) * timestep

#correlation time
tcorr = np.array(range(int(args.total_time/timestep)+1)) * timestep # in ps
tcorr_len = len(tcorr)

JxJx = autocorr(Jx)
print('Jx done!')
JyJy = autocorr(Jy)
print('Jy done!')
JzJz = autocorr(Jz)
print('Jz done!')

cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)

jxyz_full = np.array([JxJx, JyJy, JzJz]) * metal2SIps
kxyz_full = np.array([cumsum_JxJx, cumsum_JyJy, cumsum_JzJz])

jxyz = jxyz_full[:,:tcorr_len]
kxyz = kxyz_full[:,:tcorr_len]

j = np.mean(jxyz, axis=0)
k = np.mean(kxyz, axis=0)

k_avg = np.mean(k[(tcorr > args.avg1) & (tcorr < args.avg2)]) # average kappa between 0.1 and 10 ps

fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)
fig.subplots_adjust(hspace=0.05)

ax[0].plot(tcorr, j, c='k', label=f'Average ({t[-1] / 1000:.1f} ns)')
ax[0].plot(tcorr, jxyz[0], c='C1', label='x')
ax[0].plot(tcorr, jxyz[1], c='C2', label='y')
ax[0].plot(tcorr, jxyz[2], c='C0', label='z')
ax[0].set_ylabel('C(t)' + ' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1} \mathrm{ps}^{-1})$')
ax[0].plot(tcorr, np.ones(tcorr.shape)*0,'k--')

ax[1].plot(tcorr, k, c='k')
ax[1].plot(tcorr, kxyz[0], c='C1')
ax[1].plot(tcorr, kxyz[1], c='C2')
ax[1].plot(tcorr, kxyz[2], c='C0')
ax[1].axhline(k_avg, c='k', ls='--', alpha=0.7, label=f'{k_avg:.2f} ({args.avg1:.1f} - {args.avg2:.1f} ps)')

ax[1].set_xlabel('t (ps)')
ax[1].set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')

for a in ax:
    a.set_xscale('log')
    a.legend(fancybox=False, edgecolor='black')
    a.minorticks_on()

plt.savefig(args.outfile + '.jpg', dpi=300, bbox_inches='tight')

with open("Press.txt", "r") as file:
    P = float(file.read())

if args.save_pkl:
    print('Save data in pickle format ...')
    with open(args.outfile + '.pkl', 'wb') as file:
        print('Total time for autocorrelation:', t[-1])
        print('Total time for plotting:', tcorr[-1])
        data = [T, P, tcorr, jxyz, kxyz]
        pkl.dump(data, file)
else:
    print('Save data in txt format ...')
    np.savetxt(args.outfile + '.txt', np.array([tcorr, jxyz[0], jxyz[1], jxyz[2], kxyz[0], kxyz[1], kxyz[2]]).T, header=f't (ps) Jx Jy Jz kx ky kz; T = {T:.2f} K, P = {P:.2f} GPa', fmt='%.8f')
