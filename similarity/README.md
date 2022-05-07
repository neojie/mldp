# Similarity

Post-processing scripts for the soap similiarity analysis



```
mkdir soap
cp in.lammps_360pv108w /u/scratch/j/jd848/360pv108w_75g/2500/soap
lmp -in in.lammps_360pv108w  # better sumbit job
python ~/script/mldp/similarity/quick_merge.py
```



`quick_merge.py` merges the soap similarity analyzed file of each element into one .xyz file => no speicial dependencies



`stat2_pv_20g.py` presents an example of using Gibbs dividing surface to analyze the interface, routines in the `stat_lib.py` are called.  => speicial dependencies : MDAanalysis, lmfit, ase 

`stat_lib_mass.py`  collection of scripts, mass as the weighting factor





==TODO==

1- merge stat_lib.py and stat_lib_mass.py

weight can be specify by the user

2- any amounts elements 

3- flexibility of modulate the 





mass as weight 

bash_sub_mass.py

stat_lib_mass.py

stat_mass.py
