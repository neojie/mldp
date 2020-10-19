#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 09:55:29 2020

collection of functions used

@author: jiedeng
"""
import os
import re
import shutil
import glob

def get_etot_outcar(outcar):
    etot = []
    fo      = open(outcar)
    for line in fo: # https://www.codecademy.com/forum_questions/51777a4d81ce6f1632002c30, file object is iterable
        if 'free  energy   TOTEN' in line:
#            print(line)
            energy = float(line.split()[-2])
#            print(energy)
            etot.append(energy)
    return etot

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


def load_paths(inputpath,level='vasp',output=False):
    """
    load paths and modify path to level required
    
    inputpath: input path, a file name storing input paths, str
    level: level could be vasp, recal, deepmd
    output: write output as a single file
    
    return
    ------
    paths: list of str
    
    Note
    -----
    path should be:
    xx/recal/deepmd*, xx/recal, xx
    or
    xx/recal/deepmd*  relax_step, relax_paths
    
    xx is vasp level
    deepmd* could be deepmd, deepmd_relax, deepmd_relax2, deepmd_relax, etc.
    

    """
    fp          = open(inputpath,'r')
    folders_org = fp.readlines()
    paths     = []
    relax_paths = []
    relax_step  = []
    fp.close()
    relax = False
    for i in range(len(folders_org)):        
        #skip #xxx and empty
        if folders_org[i].split() == [] or  folders_org[i][0] == '#' or folders_org[i][:1] == '\n':
            pass
        else:
            # remove \n in the end
            tmp = folders_org[i].replace('\n','')
            # remove # 
            tmp = tmp.split('#')[0]
            # remove " and ,                     
            tmp= tmp.replace('"','').replace(',','').split()
            
            if len(tmp) == 3:
                relax = True
                tmp2 = _load_paths_helper(tmp[0],level)
                if not(tmp2 in paths):
                    paths.append(tmp2)
                    
                if tmp[2] == 'self':
                    relax_paths.append(paths[-1])
                else:
                    relax_paths.append(tmp[2])
                relax_step.append(int(tmp[1]))
                
            elif len(tmp) == 1:
                tmp2 = _load_paths_helper(tmp[0],level)
                if not(tmp2 in paths):
                    paths.append(tmp2)
            else:
                print('????????',tmp,'length is not 3 or 1, skip!')
                
    
    if relax_paths != []:
        paths = [paths,relax_step, relax_paths]
    
    if output:
        from datetime import date
        today = date.today()
        out = open('paths-{0}-{1}-{2}'.format(today.year,today.month,today.day),'w')
        if relax is True:
            for path in paths[0]:
                out.write(path)
                out.write('\n')
        else:
            for path in paths:
                out.write(path)
                out.write('\n')
        out.close()            
    return paths

def _load_paths_helper(path,level):
    """
    if '/recal' exist
    use recal to divide the path
    if recal does not exist
    must be vasp level
    """
    if '/recal' in path:
        vasp_path = path.split('/recal')[0]
    else:
        vasp_path = path
    
    if level == 'vasp':
        return vasp_path
    if level == 'recal':
        return os.path.join(vasp_path,'recal')
    if 'deepmd' in level:
        return os.path.join(os.path.join(vasp_path,'recal'),level)
        



########merge outcar######
def header(subfolders):
    """
    generate header of outcar
    """
    str1 = 'nsw_sel= '
    nsw_sel    = [os.path.basename(os.path.normpath(subfolder)) for subfolder in subfolders]
    nsw_sel_int  = [int(i) for i in nsw_sel]; nsw_sel_int.sort()
    nsw_sel      = [str(i) for i in nsw_sel_int]
    for i in nsw_sel:
        str1 = str1 + i + ' '
    str2 = 'nsw_tot= ' + str(len(subfolders))       
    string = '#'*100 + '\n' + str1 + '\n' + str2 + '\n' + '#'*100
    return string, nsw_sel


def merge(path):
    outcar     = os.path.join(path,'OUTCAR')
    if os.path.exists(outcar):
        print("OUTCAR exists, skip")
    else:
        print("No OUTCAR, building OUTCAR")
        subfolders = [f.path for f in os.scandir(path) 
                      if f.is_dir() and os.path.exists(os.path.join(f.path,'OUTCAR'))]   # for future usage, only depends on  outcarXXX
        string, nsw_sel     = header(subfolders)
        os.system("echo '{0}' >> {1}".format(string, outcar))
        for nsw_i in nsw_sel:
            outcarInfolder = os.path.join(os.path.join(path,str(nsw_i)),'OUTCAR')    
            os.system("echo '%nsw_i= {0}' >> {1}".format(nsw_i, outcar))
            os.system("cat {0} >> {1}".format(outcarInfolder, outcar)) 
    
def merge_sel(path,OUTCAR):
    """
    path : str, path contains subfoders of OUTCAR runs
    OUTCAR : str, name of merged OUTCAR output
    
    selectively merge outcar
    removing the redundant content in OUTCAR
    """
    outcar     = os.path.join(path,OUTCAR)
    
    if os.path.exists(outcar):
        print("OUTCAR_rleax exists, skip")
    else:
        print("No OUTCAR, building {0}".format(OUTCAR))
        subfolders = [f.path for f in os.scandir(path) 
                      if f.is_dir() and os.path.exists(os.path.join(f.path,'OUTCAR'))]   # for future usage, only depends on  outcarXXX
        string, nsw_sel  = header(subfolders)
        
        os.system("echo '{0}' >> {1}".format(string, outcar))
        
        if len(subfolders) > 1 :
            line_num = 0
            with open(os.path.join(subfolders[0],'OUTCAR')) as f:
                for line in f:
                    if 'Iteration' in line:
                        break
                    else:
                        line_num += 1
        for nsw_i in nsw_sel:
            outcarInfolder = os.path.join(os.path.join(path,str(nsw_i)),'OUTCAR')    
            os.system("echo '%nsw_i= {0}' >> {1}".format(nsw_i, outcar))
            if nsw_i  == nsw_sel[0]:
                os.system("cat {0} >> {1}".format(outcarInfolder, outcar)) 
            else:
                os.system("tail --line=+{0} {1} >> {2}".format(line_num-1, outcarInfolder, outcar)) 


def remove_recal_traj_files(path):
    subfolders = [f.path for f in os.scandir(path) if f.is_dir() and os.path.exists(os.path.join(f.path,'OUTCAR'))] 
    for subfolder in subfolders:
        shutil.rmtree(subfolder)
    print("remove XDATCAR, *out*, dsq*")
    wild_cards = ['XDATCAR','*out*','dsq*', 'job_*','*sh','*tsv']
    for wild_card in wild_cards:
        target =  os.path.join(path,wild_card)
        for file_to_remove in glob.glob(target):
            os.remove(file_to_remove)
