#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

TEMP       = 10000 # Kelvin
NUM_ATOMS  = 160 
ATOMS_MASS = 16   # in g/mol

# Constants
kB          = 1.3806488e-23 # J/K
H_MASS      = 1.673534e-27 # Mass of hydrogen atom (Kg)
ATOMS_MASS  = ATOMS_MASS*H_MASS # Kg
m_per_sec2A_per_fs = 1e-5 

# set-up
PEAK_SPEED  = np.sqrt(2 * kB *TEMP/ATOMS_MASS) # peak velocity (m/s), derived from derivative of f(v) , Velocity is a vector, speed (m/s) is magnitude of that vector. 
MAX_SPEED   = round(PEAK_SPEED*3.5) # Max velocity we set 


def MaxwellSpeedDist(s, m, T):
    # See wikipedia: https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    s2 = s*s # s speed
    twokT = 2 * kB * T
    coeff = 4 * np.pi * (m / (np.pi * twokT)) ** (3./2.)
    return coeff * s2 * np.exp(-m * s2 / (twokT))

    
def theorySpeedDist():
    tsd = list()  # probability
    spds = list() # velocity distribution
    for s in range(0, MAX_SPEED, 1):
        tsd.append(MaxwellSpeedDist(s, ATOMS_MASS, TEMP))
        spds.append(s)
    
    pctCovered = sum(tsd)
    print("Sanity check: ", pctCovered, pctCovered >= 0.99)
    return spds, tsd
    
def MCSpeedDist(maxProb):
    count = 0
    mcsd = list()
    while (len(mcsd) < NUM_ATOMS):
        s = np.random.random() * MAX_SPEED
        prob = MaxwellSpeedDist(s, ATOMS_MASS, TEMP)
        
        if np.random.random() * maxProb <= prob:
            mcsd.append(s)
        count += 1
    
    print("MC Loops: %d, Loops/Particle: %f" % (count, float(count)/NUM_ATOMS))
    return mcsd
    
spds, tsd = theorySpeedDist()
mcsd      = MCSpeedDist(max(tsd))

plt.figure()
plt.plot(np.array(spds)*m_per_sec2A_per_fs,tsd)
plt.xlabel('speed magnitude (A/fs)')
plt.ylabel('probabiltiy')

plt.figure()
plt.plot(np.array(mcsd)*m_per_sec2A_per_fs)
plt.ylabel('speed (A/fs)')
plt.xlabel('atom id')
