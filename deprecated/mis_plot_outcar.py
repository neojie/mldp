import numpy as np

from shared_functions import get_etot_outcar
import os
import matplotlib.pyplot as plt

outcar = os.path.join('.','OUTCAR')
#outcar = '/Users/jiedeng/Documents/tmp/jd848/project_folder/ml/si-16/cont2/r3/OUTCAR'
etot = get_etot_outcar(outcar)
plt.figure()
plt.plot(etot,'o',markersize = .5)
plt.xlabel('step')
plt.ylabel('TOTEN (eV)')
plt.show()
