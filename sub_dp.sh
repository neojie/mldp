#!/bin/bash

#SBATCH --job-name=deep_learn
#SBATCH --output=dp_train.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:k80:2
#SBATCH --partition=gpu
#SBATCH --time=01:30:00

module restore cuda90
source activate tf210
cd	/gpfs/loomis/project/kklee/jd848/pv+hf/5k/60g/r3-80atoms/recal/							;	python ~/script/dpkit/stat_v2.py -de deepmd_scale -m /gpfs/loomis/project/kklee/jd848/pv+hf/dp-train/scale2/dif_p5/pv-gpu.pb -e 0.5 -d all
cd 	/gpfs/loomis/project/kklee/jd848/pv+hf/5k/60g/r3-80atoms/recal/							;	python ~/script/dpkit/stat_v2.py -de deepmd_scale -m /gpfs/loomis/project/kklee/jd848/pv+hf/dp-train/scale2/dif_p5/pv-gpu.pb -e 0.5 -d all
cd /home/jd848/script/dpkit									
python stat_v2.py -ip folders/all_deepmd_scale -e 0.5 -m /gpfs/loomis/project/kklee/jd848/pv+hf/dp-train/scale2/dif_p5/pv-gpu.pb -de deepmd_scale									
