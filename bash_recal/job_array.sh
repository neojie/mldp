#!/bin/bash
#$ -cwd
#$ -o job.$JOB_ID
#$ -j y
#$ -pe shared 1
#$ -l h_rt=24:00:00
#  Email address to notify
#$ -M $USER@mail
#  Notify at beginning and end of job
#$ -m bea
source /u/local/apps/anaconda3/etc/profile.d/conda.sh
conda activate tf2gpu
python ~/script/mldp/job_array.py -mj 2290 > job.log
