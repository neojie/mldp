#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 18:05:07 2022

@author: jiedeng
"""


#import os
from subprocess import call
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--range',"-r", type=int,nargs="+",help="range of index starting from 0")
parser.add_argument('--vasp_cmd','-v',type=str,default="qsub ~/script/sub/vasp/sub_vasp_std_1cpu.sh",help="vasp command")
args   = parser.parse_args()


call("qstat -u jd848 > /u/home/j/jd848/script/jobs",shell=True)

fp = open('/u/home/j/jd848/script/jobs')
jobs=fp.readlines()
fp.close()

if os.path.exists('tmp'):
    call("rm tmp",shell=True)



beg, end = min(args.range), max(args.range)
    
try:
    job = jobs[2:][0]
    for job in jobs[2:]:
        jobid  = job.split()[0]
        call("qstat -j {0} | grep cwd ".format(jobid) + "| awk '{print $2} '>>tmp",shell=True)
except:
    print("!!No jobs running!!")
    
print("running job paths stored in tmp")

cwd=os.path.abspath(os.path.curdir)

dirs = list(range(beg,end))
refined_dirs = []

for i in dirs:
    # check if i in the tmp file
    
    path = os.path.join(cwd, str(i))
    path_ = path+'\n'
    with open('tmp') as running_job:
        if path_ in running_job.read():
            print("{0} is running".format(i))
        else:
            refined_dirs.append(path)

print("the following directories is current not running")
print(refined_dirs)

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
    for i in range(int(5)):   #somehow for some case there is 'Voluntary context switches' in between?!
        line = next(txt)
        if 'Voluntary context switches' in line:
#            print('OUTCAR is done')
            return True
    return False

                
for path in refined_dirs:
    
    print(path)
    outcar = os.path.join(path,'OUTCAR')
    flag = False
    if os.path.exists(outcar):
        try:
            outcar_done = check_outcar_done(outcar)
        except:
            outcar_done = check_outcar_done_slow(outcar)
        if outcar_done:
            flag = True 
    if not flag:
        call("cd {0};{1}".format(path,args.vasp_cmd),shell=True)
    else:
        print("VASP calculcation already done")