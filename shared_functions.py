#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 09:55:29 2020

collection of functions used

@author: jiedeng
"""
import os
#import numpy as np
import re

def check_outcar_done_slow(outcar):
    """
    slowly check if outcar is done
    require large memory if OUTCAR is large
    """
    fp  =  open(outcar)
    lines =fp.readlines() 

    string = "Voluntary context switches"
    if string in lines[-1]:
        fp.close()
        return True
    fp.close()
    return False
    

def check_outcar_done(outcar):
    """
    quickly check if outcar is done
    sometimes does not work
    """
    txt  = reverse_readline(outcar)
    for i in range(int(5)):
        line = next(txt)
        if 'Voluntary context switches' in line:
#            print('OUTCAR is done')
            return True
    return False

def check_incar(incar):
    """
    check if INCAR temperature setting correct
    return the correct incar file
    
    """
    kb = 8.617333262e-5
    with open(incar) as incar_file:
        for line in incar_file:
            if 'SIGMA' in line:
               sigma = float(line.split('SIGMA')[1].split('=')[1].split()[0])
            if 'TEBEG' in line:
               tebeg = float(re.sub('[^0-9]', '', line.split('=')[1]))
    if abs(kb*tebeg - sigma)> 0.01:
        print(incar,' SIGMA wrong')
        raise ValueError()
    else:
        sel_incar='INCAR_'+str(int(round(tebeg/1000)))+'k'
        return sel_incar
    

def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order
    read big file in reverse order is not trivial, library avaialable but not a optimum choice
    source: https://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python
    """

    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


def load_paths(inputpath):
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    fp.close()
    for i in range(len(folders_org)):
        if '#' in folders_org[i] or folders_org[i] == '\n':
            pass
        
        else:
            tmp = folders_org[i].replace('\n','')
            paths.append(tmp.replace('"','').replace(',','').replace('deepmd',''))
    return paths


