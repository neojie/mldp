
"""
batch submit jobs
"""
from pymatgen.io.vasp import Poscar
import numpy as np
import os
import glob
from shutil import copyfile
import argparse

# check if exist 

parser = argparse.ArgumentParser()
parser.add_argument("--vols","-v",help="target volume, in form of xx-xx-xx")
parser.add_argument("--prefix","-p",help="prefix, in form of r1, r3, r1 usually denotes solid and r3 liquid")

args   = parser.parse_args()

###
vols = np.array([float(i) for i in args.vols.split('-')])
#vols = np.array([1500, 1114.5, 939.1, 738.63, 518.76])/2
print("target volumes are:",vols)
###

def poscar_scale (poscar_in, poscar_out, scale) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    if 'D' == lines[7][0] or 'd' == lines[7][0]: 
        lines = poscar_scale_direct(lines, scale)
    elif 'C' == lines[7][0] or 'c' == lines[7][0] :
        lines = poscar_scale_cartesian(lines, scale)
    else :
        raise RuntimeError("Unknow poscar style at line 7: %s" % lines[7])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))

def poscar_scale_direct (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    pscale = float(lines[1])
    pscale = pscale * scale
    lines[1] = str(pscale) + "\n"
    return lines

def poscar_scale_cartesian (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    # scale box
    for ii in range(2,5) :
        boxl = lines[ii].split()
        boxv = [float(ii) for ii in boxl]
        boxv = np.array(boxv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (boxv[0], boxv[1], boxv[2])
    # scale coord
    for ii in range(8, 8+numb_atoms) :
        cl = lines[ii].split()
        cv = [float(ii) for ii in cl]
        cv = np.array(cv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (cv[0], cv[1], cv[2])
    return lines    

def poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def build_folder(vol):
    scale = (vol/vol0)**(1/3)
    folder_name = os.path.join(init_path,args.prefix + '-'+str(int(vol)))
    
    if os.path.isdir(folder_name) and os.path.isdir(os.path.join(folder_name,'POSCAR')):
        print(folder_name, 'already exist')
        pass
    else:
        os.mkdir(folder_name)
        poscar_out = os.path.join(folder_name,'POSCAR')
        poscar_scale(os.path.join(init_path,'POSCAR'), poscar_out,scale)
        copyfile(os.path.join(init_path,'POTCAR'),os.path.join(folder_name,'POTCAR'))
        copyfile(os.path.join(init_path,'INCAR'),os.path.join(folder_name,'INCAR'))
        copyfile(os.path.join(init_path,'KPOINTS'),os.path.join(folder_name,'KPOINTS'))
        copyfile(os.path.join(init_path,'sub_job.sh'),os.path.join(folder_name,'sub_job.sh'))
    return folder_name

init_path = os.getcwd()
poscar_in = os.path.join(init_path,'POSCAR')
structure = Poscar.from_file(poscar_in).structure
vol0      = structure.volume

existing_runs = glob.glob('r*')
existing_vols = []
commands = []

# folders starting from r
for folder in existing_runs:
    if 'POSCAR' in os.listdir(folder): 
        poscar_in_existed = Poscar.from_file(os.path.join(folder,'POSCAR'))
        structure         = poscar_in_existed.structure
        vol_existed       = structure.volume
        existing_vols.append(vol_existed)


if len(existing_vols) >0:  
    for vol in vols:
        folder = build_folder(vol)
        command = 'cd ' + folder + ';' 'sbatch -J {0} sub_job.sh'.format(vol)+ ';' + 'cd ' + init_path;
        commands.append(command)
else:
    existing_vols = np.array(existing_vols)
    for vol in vols:
        if (abs(vol - existing_vols)<100).any(): # The
            print(vol,'or closed ones existed, skipped')
        else:
            folder = build_folder(vol)
            command = 'cd ' + folder + ';' 'sbatch -J {0} sub_job.sh'.format(vol)+ ';' + 'cd ' +  init_path;
            commands.append(command)

with open(os.path.join(init_path,'commands.sh'), 'w') as fp :
    for com in commands:
        fp.write(com)
        fp.write('\n')
        


