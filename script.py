from dpdata import System,LabeledSystem,MultiSystems
s=System('POSCAR',fmt='poscar')
print(s)
ls=LabeledSystem('OUTCAR',fmt='outcar')
print(ls)
ls.to_deepmd_raw('deepmd')
ls.to_deepmd_npy('deepmd',set_size=1000)
