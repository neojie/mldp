#!/bin/bash
#$ -cwd
#$ -o out.$JOB_ID
#$ -j y
#$ -pe shared 1
#$ -l h_rt=24:00:00,rh7
#source /u/local/apps/anaconda3/etc/profile.d/conda.sh
. "/u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh"
#conda activate /u/home/j/jd848/project-lstixrud/dpkit2
conda activate /u/home/j/jd848/project-lstixrud/dpkit2;
. /u/local/Modules/default/init/modules.sh
module load gcc/8.3.0;module load intel/2020.4;module load cmake

export PLUMED_KERNEL=/u/home/j/jd848/project-lstixrud/plumed/lib/libplumedKernel.so
export PATH=/u/home/j/jd848/project-lstixrud/plumed/bin:$PATH
export LD_LIBRARY_PATH=/u/home/j/jd848/project-lstixrud/plumed/lib:$LD_LIBRARY_PATH
#lmp -in in.lammps
python ~/script/mldp/similiarity/stat.py
#mpirun -np 2 lmp -in in.lammps
#/u/local/compilers/intel/18.0.4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpirun -np 4 lmp -in /u/home/j/jd848/project-lstixrud/test/md/in.lammps.compressed