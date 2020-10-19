#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:04:55 2020

@author: jiedeng
"""


                         

import numpy as np

import os
import glob



def extract_sigma(line):
    tmp  =  line.split('SIGMA')[1].split('=')[1].split()[0]
    return float(tmp)



#systems = ["/gpfs/loomis/project/kklee/jd848/pv+hf/5k/60g/r3/cont3/deepmd/",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/60g/r1/cont2/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/100g/r1/cont2/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/100g/r3/cont3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/140g/r1/cont2/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/140g/r3/cont3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/300g/r3-0.8/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/300g/r3-0.9/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/300g/r3-0.95/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/5k/300g/r3-0.975/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid1/r1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid1/r3-3k/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid2/r1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid2/r4-3k-re/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid3/r1/cont2/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid3/r3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid4/r1/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/solid4/r3/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/300g/r3-0.8/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/300g/r3-0.9/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/300g/r3-0.95/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/3k/300g/r3-0.975/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/25g/r3-4k/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/60g/r3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/100g/r3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/140g/r3-4k/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/25g/r3-4k/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/60g/r3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/100g/r3/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/140g/r3-4k/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/25g/r1-0.99/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/60g/r1/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/100g/r1/cont1/deepmd",
#"/gpfs/loomis/project/kklee/jd848/pv+hf/4k/140g/r1-0.78/re2/deepmd"]
                         
systems = ["/gpfs/loomis/project/kklee/jd848/pv+hf/4k/300g/r3-0.9/deepmd",
                        "/gpfs/loomis/project/kklee/jd848/pv+hf/4k/300g/r3-0.975/deepmd",
                        "/gpfs/loomis/project/kklee/jd848/pv+hf/4k/300g/r3-0.95/deepmd",
                        "/gpfs/loomis/project/kklee/jd848/pv+hf/4k/300g/r3-0.8/deepmd"]
for path in systems:
    pwd = os.getcwd()
    os.chdir(path)
    print(path)
    #exract sigma info
    incar = path+'/../INCAR'
    with open(incar) as incar_file:
        for line in incar_file:
            if 'SIGMA' in line:
               sigma = extract_sigma(line)
               print(sigma)
    sets = glob.glob('*set*')

#    np.savetxt( 'fparam.raw', all_te)
    for seti in sets:
        print(seti)
        energy=np.load(seti+'/energy.npy')
        size = energy.size
        print("the size of the file is: ", size )
        all_te = np.ones(size)*sigma
        np.save( seti+'/fparam.npy', all_te)
    os.chdir(pwd)
print('done')
    
    





