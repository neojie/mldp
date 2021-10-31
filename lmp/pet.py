#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 08:19:41 2021


# plot pressure, temperature, and energy without entropy (default)
# output the number of atoms

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

# get nsw_tot
parser = argparse.ArgumentParser(description="extract P E T V from OUTCAR")

parser.add_argument("--beg","-b", type=int,default=0, help="begin from index")
parser.add_argument("--end","-e", type=int,default=-1, help="end at index")
args = parser.parse_args()

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

#		print "<x> = {0:f} +/- {1:f}\n".format(blockMean[-1], np.sqrt(blockVar[-1]))

		plt.tight_layout()
		plt.show()

	return v, np.sqrt(blockVar), blockMean

nsw_tot = get_nsw_tot()
target_dir = '.'
vasp_dir   = '.'
outcar = os.path.join(vasp_dir,'OUTCAR')
atoms = []
dummy = ase.io.vasp.read_vasp_out(outcar)
natoms = dummy.get_global_number_of_atoms()

print("System:", dummy.get_chemical_formula(),'\t', dummy.get_global_number_of_atoms()) # can also be done via grep NION


#property_list = []
#get_pos, get_energy = False, False
#nsw_i = 0
#ewithout = []

#if 'pet.dat' in os.listdir():
#    print("**pet.dat already exists, skip extracting P E T data**")
#    pass
#else:
#    print("**Extracting P E T data**")
#    from subprocess import call
#    call("grep 'total pressure' OUTCAR | awk '{print $4}' > p.dat",shell=True)
#    call("grep 'energy  without entropy' OUTCAR | awk '{print $4}' >e.dat",shell=True)
#    call("grep '(temperature' OUTCAR | awk '{print $6}' > t.dat",shell=True)
#    call("paste e.dat t.dat > temp",shell=True)
#    call("paste p.dat temp >pet.dat",shell=True)
#    call("rm e.dat t.dat p.dat",shell=True)

print("**Extracting P E T data**")
from subprocess import call
call("grep 'total pressure' OUTCAR | awk '{print $4}' > p.dat",shell=True)
call("grep 'energy  without entropy' OUTCAR | awk '{print $4}' >e.dat",shell=True)
call("grep '(temperature' OUTCAR | awk '{print $6}' > t.dat",shell=True)
call("paste e.dat t.dat > temp",shell=True)
call("paste p.dat temp >pet.dat",shell=True)
call("rm e.dat t.dat p.dat",shell=True)

try:
    pet= np.loadtxt('pet.dat')
except:
    ## may be text exist, one corner case that is not considered yet is
    ## each column may be of different length
    tmp=np.genfromtxt('pet.dat')
    pet=tmp[~np.isnan(tmp).any(axis=1)]


pet_be = pet[args.beg:args.end,:]
ave=np.mean(pet_be,axis=0)
print("vanilla mean of P E T V")
print(f'{ave[0]/10}\t{ave[1]}\t{ave[2]}\t{dummy.get_volume()}')

#
p, e, t = pet_be[:,0], pet_be[:,1], pet_be[:,2]

vp, bpvar, bpmean=blockAverage(p)
ve, bevar, bemean=blockAverage(e)
vt, btvar, btmean=blockAverage(t)

pmean  = np.round(p.mean(),2); pst  = np.round(p.std(),2);
emean  = np.round(e.mean(),2); est  = np.round(e.std(),2);
tmean  = np.round(t.mean(),2); tst  = np.round(t.std(),2);

# PET order
print("**P(GPA)/E(eV)/T(K) value and SR, blockaverage**")
print("%.2f \t %.2f \t %.2f" % (bpmean[-1]/10,bemean[-1],btmean[-1]))
print("%.2f \t %.2f \t %.2f" % (bpvar[-1]/10,bevar[-1],btvar[-1]))


import matplotlib.pyplot as plt
fig,ax = plt.subplots(3,1,figsize=(5,10),sharex=True)
ax[0].plot(pet[:,0]/10,label='P (GPa)')
ax[0].plot(args.beg, pet[args.beg,0]/10, 'ko')
if args.end<0:
    end = args.end + nsw_tot
else:
    end = args.end
ax[0].plot(end, pet[args.end,0]/10, 'ko')
ax[1].plot(pet[:,1],label='E (eV)')
ax[2].plot(pet[:,2],label='T (K)')
ax[2].set_xlabel("Step")
ax[0].legend();ax[0].grid()
ax[1].legend();ax[1].grid()
ax[2].legend();ax[2].grid()
plt.show()