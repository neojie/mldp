#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 12:18:50 2020

modules
ase    | https://wiki.fysik.dtu.dk/ase/ase/io/io.html
supercell => vasp, not good, poscar show mg si o mg si o mg si o, elements not grouped


@author: jiedeng
"""
print("--"*20)
print("supported format \nformat |description| capabilities|")
print("lammps-data | LAMMPS data file | RW")
print("vasp | VASP POSCAR/CONTCAR | RW")
print("pdb | VASP POSCAR/CONTCAR | W => ase un-supported, for ipi usage")
print("xyz | VASP POSCAR/CONTCAR | RW")
print("--"*20)
import ase
import ase.io
import ase.io.lammpsdata
import os
import numpy as np

class FortmatConverter(object):
    def __init__(self,
             infile,
             outfile,
             out_format='lmp',
             in_format=None,
             elements = None):
        if in_format == None:
            in_format = self.parse_fmt(infile)
        if in_format == 'lmp': # lmp does not contain element info
            Zs=ase.symbols.symbols2numbers(elements)
            Zs.insert(0,0)
            self.dat =ase.io.lammpsdata.read_lammps_data(
                    infile,style='atomic',Z_of_type=Zs)
        else:
            self.dat = ase.io.read(infile,format=in_format)
        self.out_format = out_format
        self.outfile = outfile

    def scale_v(self,newv):
        """
        scale based on volume
        """
        scale = (newv/self.dat.get_volume())**(1/3)
        self.scale_f(scale)


    def scale_f(self,scale):
        """
        scale based on volume
        """
        newcell = self.dat.cell*scale
        self.dat.set_cell(newcell,scale_atoms=True)

    def supercell(self,multiplier):
        """
        multiplier: list like, scale in x,y,z
        ----
        https://notebook.community/bkoz37/labutil/samples/lab1_samples/Ase-structure
        """
        M = np.zeros((3,3))
        M[0][0], M[1][1], M[2][2] = multiplier[0], multiplier[1], multiplier[2]
        from ase.build import make_supercell
        self.super_dat = make_supercell(self.dat, M)

    def output(self,data2output,sort):
        if self.out_format == 'lmp':
            ## workaround because ase.io.write('tmp',pos,format = 'lammps-data',direct=True)  # ase output results disorded!!!
            ##
            ase.io.write('tmp_used_by_format_converter',data2output,format = 'vasp',direct=True,vasp5=True)
            import dpdata
            ls=dpdata.System('tmp_used_by_format_converter',fmt='vasp/poscar')
            ls.to_lammps_lmp(self.outfile)
            os.remove('tmp_used_by_format_converter')
        elif self.out_format == 'vasp':
           ## for supercell, if not sort, elements not grouped, workaround is to sort and specify the elements in ascending order
            ase.io.write(self.outfile,data2output,format = self.out_format,vasp5=True,direct=True,sort=sort)
        elif self.out_format == 'pdb':
            ase.io.write('tmp_used_by_format_converter',data2output,format = 'vasp',direct=True,vasp5=True)
            pos2pdb('tmp_used_by_format_converter',self.outfile)
            os.remove('tmp_used_by_format_converter')
        else:
            ase.io.write(self.outfile,data2output,format = self.out_format)
        print("done")

def pos2pdb(poscar,outfile):
    """
    for ipi
    """
    datain = open(poscar, 'r')

    datain.readline()                 #POSCAR comment line, useless to me.

    scale=float(datain.readline())   	   #assuming ONE value
    (a1, a2, a3) = datain.readline().split()   #read in unit cell vectors,
    (b1, b2, b3) = datain.readline().split()   #but there are strings
    (c1, c2, c3) = datain.readline().split()

    a1=float(a1); a2=float(a2); a3=float(a3);  #Now I have numbers.
    b1=float(b1); b2=float(b2); b3=float(b3);
    c1=float(c1); c2=float(c2); c3=float(c3);
    eles=datain.readline()
    symbols = eles.split()
    atoms = datain.readline().split()	   #List containing how many of each.

    tmp = datain.readline().split()		#tmp need to figure out if this is...
    					#Selective dynamics, or not.

    if (tmp[0][0] == 'S'):
        SD = True
    else:
        SD = False

    if (SD):
        CorD = datain.readline().split()
    else:
        CorD = tmp

    dataout = open(outfile, 'w')
    dataout.write('TITLE cell{angstrom} positions{angstrom}\n')
    line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a1,b2,c3,90.00,90.00,90.00,' P',1)
    #line = 'CRYST1   {0:2.3f}   {1:2.3f}   {2:2.3f}  {3:2.2f}  {4:2.2f}  {5:2.2f} P 1           1\n'.format(a1,b2,c3,90.00,90.00,90.00)
    dataout.write(line)
    count = 0
    for i in range(len(atoms)):		#loop over different atoms
        symbol = symbols[i]
        for j in range(int(atoms[i])):	#loop over different posistions
            count = count+1
            tmp = datain.readline().split()
            a=float(tmp[0]); b=float(tmp[1]); c=float(tmp[2]);

            if (CorD[0][0] == 'D'):		#if DLV calculate cart. pos.
                x = scale * (a1*a + b1*b + c1*c)
                y = scale * (a2*a + b2*b + c2*c)
                z = scale * (a3*a + b3*b + c3*c)
            else:				#are ready in cart.
                x = a
                y = b
                z = c
            # format required by /u/home/j/jd848/.conda/envs/tf2cpu/lib/python3.7/site-packages/i_PI-2.0-py3.7.egg/ipi/utils/io/backends/io_pdb.py
            S="ATOM  %5i %4s%4s%6s    %8.3f%8.3f%8.3f%6.2f%6.2f           %2i\n"  % (count, symbol,'1','1', float(x), float(y), float(z),0,0,0)
    #        S="ATOM  %5d %2s   MOL     1  %8.3f%8.3f%8.3f  1.00  0.00             %2s\n" % (count, symbol, float(x), float(y), float(z),symbol)
            dataout.write(S)

    datain.close()
    dataout.close()


if __name__ == '__main__' :
    import argparse

    parser = argparse.ArgumentParser()
    ### MUST SET###
    parser.add_argument("--input","-i",default='POSCAR',help="input file")
    parser.add_argument("--output","-o",default='conf.lmp',help="output file")
    parser.add_argument("--input_fmt","-if",default='vasp',help="input format")
    parser.add_argument("--output_fmt","-of",default='lmp',help="output format")
    parser.add_argument("--scale_volume","-s",type=float,help="new volume")
    parser.add_argument("--scale_abc","-sabc",type=float,help="scale factor, new volume = volume*(scale_factor^3)")
    parser.add_argument("--supercell","-sc",default=None,type= int, nargs='+',help="supercell e.g, 2 2 2")
    parser.add_argument("--elements","-e",default=['Mg','Si','O'], nargs='+',help="required if input is lmps format")
    parser.add_argument('--sort',"-st", default=True, action='store_false',help="sort? ")
    args   = parser.parse_args()


    fc = FortmatConverter(args.input,args.output,args.output_fmt,args.input_fmt,args.elements)
    if args.scale_volume:
        fc.scale_v(args.scale_volume)
    if args.scale_abc:
        fc.scale_f(args.scale_abc)
    fc.output(fc.dat,args.sort)
    if args.supercell:
        fc.supercell(args.supercell)
        fc.output(fc.super_dat,args.sort)

#
#pos=ase.io.vasp.read_vasp('/Users/jiedeng/GD/Computation/VESTA/VData/1100.vasp')
#
##ase.build.tools.update_cell_and_positions(pos,ase.cell.Cell(np.eye(3)*10))
#
#newa = (1342.9922*1.05)**(1/3)
#newv = 1342.9922*1.05
#scale = (newv/pos.get_volume())**(1/3)
#
##newcell = ase.cell.Cell(np.eye(3)*newa)
#
#newcell = pos.cell*scale
#pos.set_cell(newcell,scale_atoms=True)
#
#ase.io.vasp.write_vasp('tmp',pos,direct=True,vasp5=True)
#ase.io.write('tmp',pos,format = 'vasp',direct=True,vasp5=True)
#import dpdata
#ls=dpdata.System('tmp',fmt='vasp/poscar')
#ls.to_lammps_lmp('/Users/jiedeng/GD/Computation/VESTA/VData/temp.lmp')
#
####



