#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 11:02:36 2021
unit 
https://lammps.sandia.gov/doc/compute_heat_flux.html

directly analyze the log.properties file with Step, Jx, Jy, Jz 

@author: jiedeng
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 22:26:44 2021
analyze output of ave/correlate
@author: jiedeng
"""


#from util import reverse_readline
#import os 
import argparse
parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
parser.add_argument("--input_file", "-i",type=str, default="log.properties",  help="log_lmp generated file")

#parser.add_argument("--num", "-n",type=int, default=200,  help=" Nrepeat in ave/correlate Nevery Nrepeat Nfreq")
#parser.add_argument("--timestep", "-ts",type=float, default=1,  help=" timestep in fs, default 1fs")
#parser.add_argument("--scale", "-s",type=float,  help=" scale to SI unit, check the log file for this value, default 1")
parser.add_argument("--average", "-a",nargs="+",type=int, help=" step window average the thermal conductivity")
parser.add_argument("--temperature", "-t",type=float,help='temperature in K')
parser.add_argument("--volume", "-v",type=float,help='volume in A3')
#parser.add_argument("--sample_rate", "-sr",type=float,help='sample rate')
parser.add_argument("-s", "--store", default=False, action='store_true', help="Defualt:  Do not save data as outfile")
parser.add_argument("-of", "--outfile",type=str,default='log.kappa', help="out file name")

args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
import glob

def autocorr(a):
    b=np.concatenate((a,np.zeros(len(a))),axis=0)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)//2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d

if args.temperature:
    T = args.temperature
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? temperature not provided, parse from in file") 
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

        
if args.volume:
    V = args.volume
else:
    infile = 'conf.lmp'
    fp = open(infile)
    ins = fp.readlines()
    print("  ?? vol not provided, parse from in conf.lmp") 
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
timestep    = 1e-3    # in ps this is different from post_corr
sample_rate = 1       # for all mgsio3 syste, sample every single step

scale = convert/kB/T/T/V*sample_rate*timestep

metal2SIps = convert/kB/T/T/V

file = args.input_file

J=np.loadtxt(file) # heat current J is saved, but heat flux (energy x velocity) autocorrelation is saved in J0Jt, heat flux/V = heat flux
J = J*V
Jx, Jy, Jz = J[:,1], J[:,2], J[:,3]

#correlation time
dt = np.array(range(len(J)))*timestep # in ps

JxJx = autocorr(Jx)
JyJy = autocorr(Jy)
JzJz = autocorr(Jz)
JJ = (JxJx + JyJy + JzJz)/3

k11 = np.trapz(JxJx)*scale
k22 = np.trapz(JyJy)*scale
k33 = np.trapz(JzJz)*scale

kappa = (k11+k22+k33)/3.0
print(' Last step kappa is {0} (W m-1 K-1): '.format( kappa))
#print('kappa is 'kappa)

import scipy.integrate
cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3

if args.store:
    header='Step (ps)    autocorr (W/m/K/ps)    kappa(W/m/K)'
    np.savetxt(args.outfile,np.concatenate(([dt],[JJ*metal2SIps],[cumsum_JJ])).T,fmt='%12.10f',header = header)

if args.average:
    window = args.average
    cumsum_JJ_sel = cumsum_JJ[range(window[0],window[1])]
#    window_mean = cumsum_JJ[range(window[0],window[1])].mean()
    print('mean kappa within {0} ps - {1} ps is {2}   {3}(W m-1 K-1): '.format(
            dt[window[0]], 
            dt[window[1]], 
            min(cumsum_JJ_sel),
            max(cumsum_JJ_sel)
            ))
    
fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
ax[0].plot(dt, JxJx*metal2SIps,label='x') # 
ax[0].plot(dt, JyJy*metal2SIps,label='y')
ax[0].plot(dt, JzJz*metal2SIps,label='z')
ax[0].plot(dt, JJ*metal2SIps,label='average')
ax[0].set_xlabel('dt (ps)')
ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
ax[0].grid(True)
ax[0].legend()
ax[0].set_xscale('log')

ax[1].plot(dt, cumsum_JxJx,label='x') # 
ax[1].plot(dt, cumsum_JyJy,label='y')
ax[1].plot(dt, cumsum_JzJz,label='z')
ax[1].plot(dt, cumsum_JJ,label='average')

ax[1].set_xlabel('dt (ps)')
ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
ax[1].grid(True)
ax[1].legend()
ax[1].set_xscale('log')
plt.show()

#plt.show()


"""
### more straighyforward way
np.sqrt(J[:,1:])

j=np.loadtxt(file)[:,1:]
j_si = j*ev2j/(ps2s*(A2m**2)) # W/m^2

jxjx = autocorr(j_si[:,0])
jyjy = autocorr(j_si[:,1])
jzjz = autocorr(j_si[:,2])

deno = 3*kB*(T**2)
B= V*(A2m**3)/deno

jj = (jxjx + jyjy + jzjz)
C = B*jj*ps2s

plt.figure()
plt.plot(C)
plt.xscale('log')
plt.grid()



fig,ax = plt.subplots(2,1,figsize=(6,10),sharex=True)
#ax[0].plot(dt, JxJx*metal2SIps,label='x') # 
#ax[0].plot(dt, JyJy*metal2SIps,label='y')
#ax[0].plot(dt, JzJz*metal2SIps,label='z')
ax[0].plot(dt, JJ*metal2SIps,label='average')
ax[0].plot(dt,C,label='C')
ax[0].set_xlabel('dt (ps)')
ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
ax[0].grid(True)
ax[0].legend()
ax[0].set_xscale('log')


trun = dt[-1]*10*5
tcorr = 2e-2


jj_bracket=(j_si**2).sum().sum()/len(j_si)  # x + y + z 

tmp = tcorr*2/trun*(jj_bracket**2)*len(j_si)
tmp**.5 *timestep*1e-12*B


(tcorr*2/trun)**.5
"""