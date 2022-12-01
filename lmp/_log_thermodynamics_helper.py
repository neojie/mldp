#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:29:41 2022

@author: jiedeng
"""
from lammps_logfile import File

def check(dat):
    for ele in dat:
        if str(ele).replace('.','').replace('e','').replace('d','').replace('+','').replace('-','').isdigit(): # data may be in int or float, may contain . , e , d , + , -
            pass
        else:
            return False
    return True

def select(dat):
    selected_idx = []
    for i  in range(len(dat)):
        if str(dat[i]).isdigit():
            selected_idx.append(i)
    return selected_idx


def _extract_data_from_log(input_file, VARS,beg_idx,run_num=-1):
    
    try:
        log = File(input_file)
    except:
        from subprocess import call
        call("cp {0} {1}".format(input_file, input_file+'tmp'), shell=True) # do not change in the original file, better for checking on running sinulation
        call("sed -i 's/style restartinfo set but has//' {0}".format(input_file+'tmp'), shell=True)
        log = File(input_file+'tmp')
        call("rm {0}".format(input_file+'tmp'), shell=True)

    # run_num  = 0
    run_num  = -1
    x   = log.get('Step',run_num=run_num)
    ys0  = [log.get(y,run_num=run_num) for y in VARS]
    
    
    Step = log.get('Step',run_num=run_num)
    if not check(Step):
    #    print('    Step col is:', Step[:10])
        print('**data messed up**')
        selected_idx = select(Step)
        
        x = (x[selected_idx]).astype(float)
        ys0 = [(y[selected_idx]).astype(float) for y in ys0]
        print('**Fixed**')
    
    # output data as a dictionary with the 
    
    dic = {}
    dic['Step'] = x[beg_idx:]
    # print(ys0)

    for ii in range(len(ys0)):
        VAR = VARS[ii]
        dic[VAR] = ys0[ii][beg_idx:]
    return dic


