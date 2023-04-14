#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:01:06 2022

nvt works well
nve not tested yet
npt does not work, give a cp value too large

@author: jiedeng
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from _log_thermodynamics_helper import _extract_data_from_log

parser = argparse.ArgumentParser(description="Calculate thermodynamic properties")
parser.add_argument("--input_file","-i", type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("-p", "--plot", default=False, action='store_true', help="Defualt: plot")
parser.add_argument("-n", "--natoms", default=160, type=int, help="natoms")
parser.add_argument("--run_num", "-r", default=-1, type=int, help="natoms")
parser.add_argument("--volume", "-v",type=float,help='volume in A3')
parser.add_argument("--ensemble", "-en",type=str,help='ensemble, nvt, nve, npt')
parser.add_argument("-b", "--beg_idx", type=int, default=0, help="Defualt: plot starting from index, sometimes, first few frames are not physical and should be excluded")

args = parser.parse_args()

# check if it is nvt nve npt
# if nvt, we need 
# extract data from the log.lammps

# get volume from starting file
if args.volume:
    volume = args.volume
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
        
    volume = (xhi - xlo)*(yhi - ylo)*(zhi - zlo)
    print(' ** volume = ', volume)


nve_vars = ['Press','Temp', 'PotEng','KinEng']
nvt_vars = ['Press','Temp', 'PotEng','KinEng']
npt_vars = ['Press','Temp', 'PotEng','Volume','Enthalpy']

ensemble = args.ensemble
if ensemble == 'nve':
    VARS = nve_vars
elif ensemble == 'nvt':
    VARS = nvt_vars
elif ensemble == 'npt':
    VARS = npt_vars



dic = _extract_data_from_log(args.input_file, VARS,args.beg_idx,run_num=args.run_num)

if (ensemble == 'nve') or (ensemble == 'nvt'):
    dic['Volume'] = volume
# constants 
kb =  8.617e-5 # eV/K
kb_natoms = kb*args.natoms
# global variable UPPER case
NATOMS = args.natoms
KB_NATOMS = kb*NATOMS
j2eV = 6.242e18
ev_A3_to_Pa = 1/j2eV/(1e-30)
kbar_Ang3_to_eV = 1e5*1e-30*j2eV


# mean variables
def fluct(data):
    tmp=data - data.mean()
    out = (tmp**2).mean()
    return out



### NPT
def npt_cp(enthalpy,temp):
    # cp
    cp = []        
    for i in range(len(enthalpy)):
        cp_i = fluct(enthalpy[i:])/(temp**2)/kb/KB_NATOMS
        cp.append(cp_i)

    return np.array(cp)
    
def npt_beta_t(volume,temp):    
    # beta_t : isothermal compressibility
    beta_t = []
    for i in range(len(volume)):
        beta_t_i = fluct(volume[i:])/temp/kb/volume.mean()
        beta_t.append(beta_t_i/ev_A3_to_Pa)
    return np.array(beta_t)
    
def npt_alpha(enthalpy,volume,temp):    
    deltaV = volume - np.mean(volume)
    deltaenthalpy = enthalpy - np.mean(enthalpy)
    
    # alpha : thermal expansitivity
    alpha = []
    for i in range(len(deltaV)):
        term1 = (deltaV[i:]*1e5*deltaenthalpy[i:]).mean()/kb/temp/temp/volume.mean()
        alpha.append(term1*ev_A3_to_Pa)
    return np.array(alpha)
    
def npt_gamma_v(alpha,beta_t):
    # gamma_v : (dP/dT)_v
    gamma_v = alpha/beta_t
    return gamma_v


### NVT
def nvt_cv(poteng,temp):
    # cp
    cv = []   
    for i in range(len(poteng)):
        cv_i = fluct(poteng[i:])/(temp**2)/kb/KB_NATOMS
        cv.append(cv_i)
    return np.array(cv)
    
def nvt_gamma_v(press,poteng,temp,volume):
    # gamma_v : (dP/dT)_v
    deltaP = press - np.mean(press)
    deltaPot = poteng - np.mean(poteng)    
    gamma_v = [] # dP/dT
    for i in range(len(deltaP)):
    
        term1 = (deltaP[i:]*1e5*deltaPot[i:]/j2eV).mean()/(kb/j2eV)/temp/temp + NATOMS/volume*kb*ev_A3_to_Pa
        gamma_v.append(term1)
    return np.array(gamma_v)

### NVE
def nve_cv(poteng,temp):
    cv = []   
    for i in range(len(poteng)):
        cv_i = 1.5/(1-fluct(poteng[i:])/1.5/kb/KB_NATOMS/temp/temp) 
        cv.append(cv_i)
    return np.array(cv)
    
def nve_gamma_v(press,poteng,temp,volume):
    # gamma_v : (dP/dT)_v
    deltaP = press - np.mean(press)
    deltaPot = poteng - np.mean(poteng) 
    cv_mean = nve_cv(poteng,temp)[0]  # 
    term2 = 2/3*cv_mean*KB_NATOMS/volume
    gamma_v = [] # dP/dT
    
    for i in range(len(deltaP)):
        
        term1 =(1-(deltaP[i:]*1e5*deltaPot[i:]).mean()*volume*1e-30*j2eV/kb/KB_NATOMS/temp/temp)    
        gamma_v.append(term1*term2*ev_A3_to_Pa)

    return np.array(gamma_v)

def plot(out):
    num=len(out)
    horizontal_scale = 3*num
    fig,ax = plt.subplots(num,1,figsize=(6,horizontal_scale),sharex=True,sharey=False)
    i = 0
    for key in out.keys():
        ax[i].plot(dic['Step'],out[key],label=key,marker='*', alpha=0.6, linestyle='-')
        ax[i].legend();ax[i].grid()
        i += 1
    plt.show()


if ensemble == 'nve':
    cv      = nvt_cv(dic['PotEng'], dic['Temp'].mean())
    gamma_v = nvt_gamma_v(dic['Press'],dic['PotEng'], dic['Temp'].mean(),dic['Volume'])
    print('cv and gamma_v')
    print(cv[0])
    print(gamma_v[0])
    out = {'cv':cv,'gamma_v':gamma_v}
    plot()
    
elif ensemble == 'nvt':
    cv      = nvt_cv(dic['PotEng'], dic['Temp'].mean())
    gamma_v = nvt_gamma_v(dic['Press'],dic['PotEng'], dic['Temp'].mean(),dic['Volume'])
    print('cv and gamma_v')
    print(cv[0])
    print(gamma_v[0])
    out = {'cv':cv,'gamma_v':gamma_v}

elif ensemble == 'npt':
    # calculate enthapy
    enthalpy = dic['PotEng'] + dic['Press']*dic['Volume']*kbar_Ang3_to_eV
    cp = npt_cp(enthalpy,dic['Temp'].mean())
    beta_t = npt_beta_t(dic['Volume'],dic['Temp'].mean())
    alpha = npt_alpha(enthalpy,dic['Volume'],dic['Temp'].mean())
    gamma_v = npt_gamma_v(alpha,beta_t)
    print('cp and gamma_v')
    print(cp[0])
    print(gamma_v[0])
    out = {'cp':cp,'gamma_v':gamma_v, 'enthalpy':enthalpy, 'beta_t':beta_t, 'alpha':alpha}
else:
    print('No implementation for this ensemble yet')


if args.plot:
    plot(out)

