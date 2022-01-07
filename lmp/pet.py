#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 08:19:41 2021


# plot pressure, temperature, and energy without entropy (default)
# output the number of atoms
# call("grep -A 4 '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' OUTCAR | grep 'energy  without entropy' | awk '{print $4}' >e.dat",shell=True)
# call("grep -A 4 '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' OUTCAR | grep 'free  energy   TOTEN' | awk '{print $5}' > etoten.dat",shell=True)

# from subprocess import call
# call("grep SCALEE OUTCAR | awk '{print $3}' > scalee.dat",shell=True)
# scalee = np.loadtxt('scalee.dat')
# if scalee < 1:
#     print("scalee = {0}< 1".format(scalee))
#     # call("grep -A 4 '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' OUTCAR | grep 'energy  without entropy' | awk '{print $4}' >e.dat",shell=True)
#     # call("grep -A 4 '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' OUTCAR | grep 'free  energy   TOTEN' | awk '{print $5}' > etoten.dat",shell=True)
#     call("grep ab-initio OSZICAR  | awk '{print $2}'>ab.dat")
#     call("grep wca OSZICAR  | awk '{print $4}'>wca.dat")

# else:
#     call("grep 'energy  without entropy' OUTCAR | awk '{print $4}' >e.dat",shell=True)
# call("grep 'total pressure' OUTCAR | awk '{print $4}' > p.dat",shell=True)
# call("grep '(temperature' OUTCAR | awk '{print $6}' > t.dat",shell=True)
# call("grep 'EENTRO' OUTCAR | awk '{print $5}' > eentro.dat",shell=True)
# call("paste e.dat t.dat > temp",shell=True)
# call("paste p.dat temp >pet.dat",shell=True)
# call("rm e.dat t.dat p.dat scalee.dat",shell=True)
            
    

# call("paste e.dat t.dat > temp",shell=True)
# call("paste p.dat temp >pet.dat",shell=True)
# call("rm e.dat t.dat p.dat scalee.dat",shell=True)
    
@author: jiedeng
"""
import os
import ase.io
import ase
import ase.io.vasp
import re
import numpy as np
import warnings
import argparse
import matplotlib.pyplot as plt

# get nsw_tot
parser = argparse.ArgumentParser(description="extract quantities from OUTCAR")

helper = 'what to show?  Default:p, ewithout,t\n'\
        'could be:\n'\
        'p : pressure\n'\
        't :  temperature\n'\
        'vol :  volume\n'\
        'toten :  free  energy   TOTEN\n'\
        'ewithout :  energy  without entropy\n'\
        'scalee :  E(ab-initio) - E(wca)'\
        'eentro :  EENTRO'
print(helper)
parser.add_argument("--beg","-b", type=int,default=0, help="begin from index")
parser.add_argument("--end","-e", type=int,default=-1, help="end at index")
parser.add_argument("-p", "--plot", default=True, action='store_false', help="Defualt: plot")
parser.add_argument("-s", "--save", default=False, action='store_true', help="Defualt: DONOT save plot")
parser.add_argument("-ep", "--extract_pet", default=True, action='store_false', help="Defualt: plot")
parser.add_argument("--list","-l",default=['p','ewithout','t'], nargs='+',help=helper)
parser.add_argument("-m", "--marker", type=str,default='', help="marker style supported by Matplotlib")
parser.add_argument("-ls", "--linestyle", type=str,default='-', help="line style supported by Matplotlib, e.g., None means no lines")
args = parser.parse_args()

scalee_command = "grep ab-initio OSZICAR | awk '{print $2}' > ab.dat ;"\
                "grep wca OSZICAR | awk '{print $4}' > wca.dat;"\
                "paste -d- ab.dat wca.dat | bc > du.dat"
# print(scalee_command)
commands = {'p':"grep 'total pressure' OUTCAR | awk '{print $4}' > p.dat",
            't' :  "grep '(temperature' OUTCAR | awk '{print $6}' > t.dat",
            'vol' : "grep 'volume of cell' OUTCAR | awk '{print $5}' > vol.dat",
            'toten' :  "grep 'free  energy   TOTEN' OUTCAR | awk '{print $5}' > toten.dat" ,
            'ewithout' : "grep 'energy  without entropy' OUTCAR | awk '{print $4}' >e.dat",
            'scalee' :  scalee_command,
            'eentro' :  "grep 'EENTRO' OUTCAR | awk '{print $5}' > eentro.dat",
            }

files = {'p':"p.dat",
        't' :  "t.dat",
        'vol' : "vol.dat",
        'toten' :  "toten.dat" ,
        'ewithout' : "e.dat",
        'scalee' :  "du.dat",
        'eentro' :  "eentro.dat",
            }

def get_nsw_tot(vasp_dir = '.'):
    if os.path.exists(os.path.join(vasp_dir,'XDATCAR')):
        xdatcar = os.path.join(vasp_dir,'XDATCAR')
        for line in reversed(list(open(xdatcar))):
            tmp = line.rstrip()
            if "Direct" in tmp:
                nsw_tot = int(re.sub("[^0-9]", "",tmp))
                break
    elif os.path.exists(os.path.join(vasp_dir,'OSZICAR')):
        warnings.warn(
            "\n ****** XDATCAR does not exist ******\n"
            + "If nsw set > 9999, we cannot extract the results \n"
            + " And you should set nsw_tot manually!"
        )
        osizcar = os.path.join(vasp_dir,'OSZICAR')
        for line in reversed(list(open(osizcar))):
            tmp = line.rstrip()
            if tmp[0] ==' ' and "T=" in tmp:
                nsw_tot = int(tmp.split()[0])
                break
    # if OSIZCAR not exist, try OUTCAR
    elif os.path.exists(os.path.join(vasp_dir,'OUTCAR')):
        warnings.warn(
            "\n ****** XDATCAR does not exist ******\n"
            + "If nsw set > 9999, we cannot extract the results \n"
            + " And you should set nsw_tot manually!"
        )
        outcar = os.path.join(vasp_dir,'OUTCAR')
        for line in reversed(list(open(outcar))):
            tmp = line.rstrip()
            if "Iteration" in tmp:
                nsw_tot = int(tmp.split()[2].split('(')[0])
                break
    elif os.path.exists(os.path.join(vasp_dir,'pet.dat')):
        try:
            pet= np.loadtxt('pet.dat')
        except:
            ## may be text exist, one corner case that is not considered yet is
            ## each column may be of different length
            tmp=np.genfromtxt('pet.dat')
            pet=tmp[~np.isnan(tmp).any(axis=1)]
        nsw_tot = len(pet)  
    return nsw_tot

def read_outcar(vasp_dir,index):
    outcar = os.path.join(vasp_dir,'OUTCAR')
    return ase.io.vasp.read_vasp_out(outcar,index)

def blockAverage(datastream, isplot=False, maxBlockSize=0):
	"""This program computes the block average of a potentially correlated timeseries "x", and
	provides error bounds for the estimated mean <x>.
	As input provide a vector or timeseries "x", and the largest block size.

	Check out writeup in the following blog posts for more:
	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty_14.html
	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html
	"""

	Nobs         = len(datastream)           # total number of observations in datastream
	minBlockSize = 1;                        # min: 1 observation/block

	if maxBlockSize == 0:
		maxBlockSize = int(Nobs/4);        # max: 4 blocs (otherwise can't calc variance)

	NumBlocks = maxBlockSize - minBlockSize   # total number of block sizes

	blockMean = np.zeros(NumBlocks)               # mean (expect to be "nearly" constant)
	blockVar  = np.zeros(NumBlocks)               # variance associated with each blockSize
	blockCtr  = 0

				#
				#  blockSize is # observations/block
				#  run them through all the possibilities
				#

	for blockSize in range(minBlockSize, maxBlockSize):

		Nblock    = int(Nobs/blockSize)               # total number of such blocks in datastream
		obsProp   = np.zeros(Nblock)                  # container for parcelling block

		# Loop to chop datastream into blocks
		# and take average
		for i in range(1,Nblock+1):

			ibeg = (i-1) * blockSize
			iend =  ibeg + blockSize
			obsProp[i-1] = np.mean(datastream[ibeg:iend])

		blockMean[blockCtr] = np.mean(obsProp)
		blockVar[blockCtr]  = np.var(obsProp)/(Nblock - 1)
		blockCtr += 1

	v = np.arange(minBlockSize,maxBlockSize)

	if isplot:

		plt.subplot(2,1,1)
		plt.plot(v, np.sqrt(blockVar),'ro-',lw=2)
		plt.xlabel('block size')
		plt.ylabel('std')

		plt.subplot(2,1,2)
		plt.errorbar(v, blockMean, np.sqrt(blockVar))
		plt.ylabel('<x>')
		plt.xlabel('block size')
		plt.tight_layout()
		plt.show()

	return v, np.sqrt(blockVar), blockMean

nsw_tot = get_nsw_tot()
    
target_dir = '.'
vasp_dir   = '.'
outcar = os.path.join(vasp_dir,'OUTCAR')
atoms = []
try:
    dummy = ase.io.vasp.read_vasp_out(outcar)
    natoms = dummy.get_global_number_of_atoms()
    formula =  dummy.get_chemical_formula()
    natoms = dummy.get_global_number_of_atoms()
    volume = dummy.get_volume()
except:
    formula= 'unknown'
    natoms = 'unknown'
    volume = 'unknown'

print("System:", formula,'\t','natoms:', natoms,'\t','nsw:', nsw_tot) # can also be done via grep NION



if args.extract_pet:
    print("**Extracting data**")
    from subprocess import call
    for li in args.list:
        call(commands[li],shell=True)
        if li == 'scalee':
            np.savetxt('du.dat',np.loadtxt('ab.dat')-np.loadtxt('wca.dat'))
else:
    print("**use existing P E T data pet.dat in the current folder**")

res = {}
means = []
vares = []
for li in args.list:
    file = files[li]
    res[li] = np.loadtxt(file)
    vp, var, ave=blockAverage(res[li][args.beg:])
    means.append(ave[-1])
    vares.append(var[-1])

# PET order
print(*args.list,sep='\t')
print(*means,sep='\t')
print(*vares,sep='\t')


def running_average(x):
    return np.flip(np.flip(x).cumsum()/(np.array(range(0,len(x)))+1))

    horizontal_scale = len(args.list)*3
    fig,ax = plt.subplots(len(args.list),1,figsize=(6,horizontal_scale),sharex=True,sharey=False)
    if len(args.list) == 1:
        ax.plot(res[args.list[0]][args.beg:],label=args.list[0],marker=args.marker, alpha=0.6,linestyle=args.linestyle)
        ax.legend();ax.grid()
    else:
        for i in range(len(args.list)):
            ax[i].plot(x[args.starting_i:], res[args.list[i]][args.beg:],label=args.list[i],marker=args.marker, alpha=0.6, linestyle=args.linestyle)
            ax[i].legend();ax[i].grid()
    plt.show()
if args.plot:
    horizontal_scale = len(args.list)*3
    fig,ax = plt.subplots(len(args.list),1,figsize=(6,horizontal_scale),sharex=True,sharey=False)
    if len(args.list) == 1:
        ax.plot(res[args.list[0]][args.beg:],label=args.list[0],marker=args.marker, alpha=0.9,linestyle=args.linestyle)
        ax.legend();ax.grid()
    else:
        for i in range(len(args.list)):
            ax[i].plot(res[args.list[i]][args.beg:],label=args.list[i],marker=args.marker, alpha=0.9, linestyle=args.linestyle)
            ax[i].legend();ax[i].grid()
    plt.show()

if args.save:
    fig.savefig('pet.png',bbox_inches='tight')
plt.show()

# if args.plot:
#     # import matplotlib.pyplot as plt
#     fig,ax = plt.subplots(3,2,figsize=(8,12),sharex=True)
#     ax[0][0].plot(pet[:,0]/10,label='P (GPa)')
#     ax[0][0].plot(pet[:,0]/10,label='P (GPa)')
#     ax[0][0].plot(args.beg, pet[args.beg,0]/10, 'ko')
#     if args.end<0:
#         end = args.end + nsw_tot
#     else:
#         end = args.end
#     ax[0][0].plot(end, pet[args.end,0]/10, 'ko')
    
#     ax[1][0].plot(pet[:,1],label='E (eV)')
#     ax[2][0].plot(pet[:,2],label='T (K)')
#     ax[2][0].set_xlabel("Step")
#     ax[0][0].legend();ax[0][0].grid()
#     ax[1][0].legend();ax[1][0].grid()
#     ax[2][0].legend();ax[2][0].grid()
    
#     ax[0][1].plot(running_average(pet[:,0]/10));ax[0][1].grid()
#     ax[1][1].plot(running_average(pet[:,1]));ax[1][1].grid()
#     ax[2][1].plot(running_average(pet[:,2]));ax[2][1].grid()
#     if args.save:
#         fig.savefig('pet.png',bbox_inches='tight')
#     plt.show()