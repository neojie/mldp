#!/bin/bash
#  Use current working directory
#$ -cwd
#$ -o ./joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#
#  Resources requested
#$ -pe dc* 1
#$ -l h_data=4096M,h_rt=24:00:00
#
#  Job is not rerunable
#$ -r n
# Initialization for mpi parallel execution
#
. "/u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh"
conda activate /u/home/j/jd848/project-lstixrud/dpkit2;
python job_array_dpdata.py -mj 3000 -d deepmd -if /u/project/ESS/lstixrud/jd848/metad/pvh/inputs/inputs_6.5k/ > job.log
