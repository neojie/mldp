# MLDP

## 1. Recalculate

## 2. Pertubation

##### Perturb systems and run simulations
#### dependencies
---------------------------
-`ase`
-`MDAnalysis`
-`dpdata`

### 2.1. Workflow
---------------------------
1. analyze rdf with `MDAnayalysis`? search for which two pairs to swap so that the short interatomic distance of the corresponding pair can be reached. `velocity.py` calculate velocity of the atoms given timestep and temperature and determine the min interatomic distance should be reached given target pressure, temperature, and timestep
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
### 2.1. Workflow for thermal conductivity calculation
1. run lammps calculation 
2. `log_lmp.py` extract the v_Jx, v_Jy, v_Jz heat current 
```python ~/script/mldp/lmp/log_lmp.py log.lammps -y v_Jx v_Jy v_Jz -s -p```
4. If ave/correlate output file is stored, which is NOT recommended for liquid since many auto-correlation needs to be done, `post_corr.py` can be used to analyze the results
5. If no ave/correlate output file, use `kappa.py` analyze the output of step 2)
```python ~/script/mldp/lmp/kappa.py -s -a 500 1500```
-a specifcy average between step 500 to 1500
 ## 4. ASAP
Fingerprint analysis using ASAP and Dscribe
#### dependencies
---------------------------
-`ASAP`
-`Dscribe`
 
## 5. Model deviation 

## 6. Scale
## 7. Util
> Written with [StackEdit](https://stackedit.io/).
