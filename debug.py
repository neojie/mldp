#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 20:04:37 2020

@author: jiedeng
"""
inputpath = 'folders_relax'
def load_paths(inputpath):
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    relax_paths = []
    relax_step  = []
    fp.close()
    for i in range(len(folders_org)):
#        print(i)
        #skip #xxx and empty
        if folders_org[i].split() == [] or '#' == folders_org[i].split()[0] or folders_org[i].split()[0] == '\n':
            pass
        else:
            # remove \n in the end
            tmp = folders_org[i].replace('\n','')
            # remove # 
            tmp = tmp.split('#')[0]
            # remove " and ,                     
            tmp= tmp.replace('"','').replace(',','').split()
            if len(tmp) == 3:
                paths.append(tmp[0].replace('deepmd',''))
                if tmp[2] == 'self':
                    relax_paths.append(paths[-1])
                else:
                    relax_paths.append(tmp[2])
                relax_step.append(int(tmp[1]))
            elif len(tmp) == 1:
                paths.append(tmp[0].replace('deepmd',''))
            else:
                print(folders_org[0],'wrong!')
    if relax_paths != []:
        paths = [paths,relax_step, relax_paths]
    return paths

paths = load_paths(inputpath)