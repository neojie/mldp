# MLDP

## 0. Preparation

create necessary directories and initial configuration

```bash
mkdir ml
cd ml
mkdir nn inputs 
mkdir 1st 2nd 3rd 4th # nth iterations
cd nn; mkdir scripts; # copy train compress test job submission scripts of both GPU and CPU versions here
cd ../..
cd inputs;
## build input files for all the endmembers, if AB, then 
mkdir A B AB
## build poscar files, stores all the possible poscars that you want to include in the training set
mkdir pos
## create a config document, where 
cd ..
vi config  
```

config should contains the paths to ml, nn, inputs, maybe specify the POTCAR configurations as well?

```bash
#some variables
#mldp path
mldp=~/script/mldp
#inputs path
inputs=/u/scratch/j/jd848/ml/inputs
#nn path
nn=/u/scratch/j/jd848/ml/nn
```

run `bash config` before proceeding

## 1. Recalculate

#### 1.1:  recalculate with VASP

1. Generate Descriptor using `ASAP`

   ```bash
   asap gen_desc -s 10 --fxyz OUTCAR soap -e -c 6 -n 6 -l 6 -g 0.44
   asap gen_desc -s 1 --fxyz npt.dump soap -e -c 6 -n 4 -l 4 -g 0.44
   ```
   
2. PCA analysis `fps` to identify frames to re-calcualte

   ```bash
   python $mldp/asap/select_frames.py -i ASAP-desc.xyz -n 70 -s 10
   python $mldp/asap/select_frames.py -i ASAP-desc.xyz -n 100
   ```
   
3. `extract_deepmd.py` with `-id` flag and index file generated in last step, consider build a `pre` folder first

   ```bash
   python $mldp/extract_deepmd.py -f OUTCAR -id index_file -st # OUTCAR contains temperature info
   python $mldp/extract_deepmd.py -f ../npt.dump -fmt dump -id ../test-frame-select-fps-n-100.index -st -t 4000
   ```
   
4. prepare a folders named `input` with `INCAR`,`KPOINTS`, `POTCAR`,`sub_vasp.sh`. Files must be tested for convergence. Also NBANDS and NELEM should be sufficient. Use `recal_dpdata.py` to recalculate selected frames

   ```bash
   python $mldp/recal_dpdata.py -d deepmd/ -if $inputs/mgofe/inputs_5000 -rv no
   python $mldp/recal_dpdata.py -d deepmd/ -if $inputs/mgofe/inputs_4000 -rv no
   ```

5. Inside`recal`folder

   `python $mldp/post_recal.py -ss $inputs/sub_vasp.sh ` 

6. Inside `recal` folder, do `python $mldp/check_nbands_nelm.py -ip all -v`

7. Inside `recal` folder, do`python $mldp/merge_out.py -o OUTCAR -r y`

8. Inside `recal` folder, remove the old deepmd folder, do 

   ````bash
   python $mldp/extract_deepmd.py -d deepmd -ttr 10000
   ````

9. `dp test`

   ```bash
   dp test -m $nn/m4/v1/pv.pb -d m4v1 -n 400
   ```

10. analyze nn and vasp

  ```bash
  python $mldp/model_dev/analysis.py -tf . -mp m4v1 -rf . -euc 10 -fuc 10 -flc 0.3
  python $mldp/model_dev/analysis.py -tf . -mp m1v2 -rf . -euc 10 -fuc 10 -flc 0.6
  python $mldp/model_dev/analysis.py -tf . -mp m2v1 -rf . -euc 10 -fuc 10 -flc 0.6 -elc 0.8
  python $mldp/model_dev/analysis.py -tf . -mp m5v1 -rf . -euc 10 -fuc 10 -flc 0.4  -elc 0.02
  ```

11. build `deepmd` based on the idx file generated and remove/keep the old deepmd

12. `dp train `



####1.2:  model deviation

1. extract frames

   ```bash
   python $mldp/extract_deepmd.py -f ./dump.0 -fmt dump -ttr 1000000 -t 3000 -st 4500
   
   ```
   
2. dp test with different models

   ```bash
   dp test -m $nn/m3/m3v1/m3v1.comp.pb -d m3v1 -n 2000
   dp test -m $nn/m1/m2/pvh4.comp.pb -d m3v2
   ```

3. analyze model deviation

   ```bash
   python $mldp/model_dev/analysis.py -tf . -mp m2-pvh4 -euc 10 
   python $mldp/model_dev/analysis.py -rf . -tf . -mp m3v1 m3v2 m3v3 m3v4 
   ```

   the upper/lower limits of force and energy RMSEs should be benchmarked with VASP runs at least ones

   idx file `idx_model_deviaton` is derived

   ```bash
   python $mldp/model_dev/analysis.py -rf . -tf . -mp m3v1 m3v2 m3v3 m3v4 -elc 0.01 -euc 2 -flc 0.4 -fuc 1
   ```

   

4. A new `dump` file that contatins a subset of the frames of the original dump file is built based on  `idx_model_deviaton`  with `dump.py`

   ```bash
   python $mldp/lmp/subset_dump.py -h
   ```

   

5. Generate Descriptor using `ASAP` on the new dump file

   ```bash
   asap
   ```

   visual inspection!

6. follow the section **1.2** for the rest of steps.

   

## 2. Pertubation

##### Perturb systems and run simulations
#### dependencies
---------------------------
-`ase`
-`MDAnalysis`
-`dpdata`

### 2.1. Workflow
---------------------------
1. analyze rdf with `MDAnayalysis` search for which two pairs to swap so that the short interatomic distance of the corresponding pair can be reached. `velocity.py` calculate velocity of the atoms given timestep and temperature and determine the min interatomic distance should be reached given target pressure, temperature, and timestep
2. simulation using good POSCAR without pertubation, set temperature and teimstep based on 1), to make the system collapse as quickly as possible
3. simulation using perturbed POSCAR
`pert.py` perturbed POSCAR , input must be vasp/poscar format, for now this file is designed for MgSiO3 only, automatically swap Si-O, Mg-O, Mg-Si
`post_pert.py` inpsect the interatomic distnace for piars in Mg-Si-O system, current cutoffs are designed for MgSiO3 up to around 1400 GPa. Support `dump`, `vasp/poscar`, `lmp` format
4. simulation with lammps
```lmp -in in.lammps``` login node usually can handle this. Do check the interatomic distance frequently. You do not want to waste time doing unnecessary runs
5. check interatomic distance.
    ```python ~/script/mldp/pert/post_pert.py -f mgsio3.dump -ft dump```
    ==WARNING: Dump file may have lost atoms. If so, corresponding frames should be deleted==
6. recal with `recal_lmp.py`
```python ~/script/mldp/recal_lmp.py -if /u/home/j/jd848/project-lstixrud/pv+hf/dp-train/lmp_run/6k/rp5/160-cpu/pert/4k_mgo_swap_p2/inputs -r 0-7```
Here step 5) output 0-6 generate interatomic distance close within the range prescribed
7. check if all vasp runs are done
```python ~/script/mldp/post_recal_lmp.py```
If not, `bash out`
8. check if all runs have sufficent nbands and NELM is good
```python ~/script/mldp/check_nbands_nelm.py -ip all```
if not, increase NBANDS, NELM in INCAR
9. Merge all vasp runs to one single `OUTCAR`
    ```python ~/script/mldp/merge_out.py -o OUTCAR -r y``` 
    ==Be cautious about the -r (remove everything) flag==
10. Build `deepmd` input file from `OUTCAR`
    ```python ~/script/mldp/extract_deepmd.py -t  -bs 1000```
    ==1000 is a random large number so that only one set is generated, -t means no test set==


## 3. Lammps
scripts used for analyzing lammps output
#### dependencies
---------------------------
-`lammps_logfile`
### 3.1. Workflow for thermal conductivity calculation
1. run lammps calculation 

2. `log_lmp.py` extract the v_Jx, v_Jy, v_Jz heat current 
    ```python ~/script/mldp/lmp/log_lmp.py -i log.lammps -y v_Jx v_Jy v_Jz -s -p```

  for multicomponent liquid system, one should also subtract the partial enthalpy term $h_a$ (Eq 4 in Deng and Stixrude, 2021) , so the correct command should be

  ```python ~/script/mldp/lmp/log_lmp.py -i log.lammps -y v_jhx v_jhy v_jhz -s -p```

  If ave/correlate output file is stored, which is NOT recommended for liquid since many auto-correlation needs to be done, `post_corr.py` can be used to analyze the results

5. If no ave/correlate output file, use `kappa.py` analyze the output of step 2)
```python ~/script/mldp/lmp/kappa.py -s -a 500 1500```
-a specifcy average between step 500 to 1500
 ## 4. ASAP
Fingerprint analysis using ASAP and Dscribe

dpkit environment works

#### dependencies
---------------------------
-`ASAP`
-`Dscribe`

## 5. Model deviation 

`extract_deepmd.py` extract configures

`dp test` with all models, `model_dev/model_dev.py ` gives an example. Do customize the code to your need

`analysis.py` analyze the `test` results

`post_model_dev.py` postprocess 



## 6. Scale
`scale`

## 7. Util
some useful routines
