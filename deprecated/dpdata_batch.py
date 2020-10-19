from dpdata import System,LabeledSystem,MultiSystems
import os


fp          = open('folders_to_merge','r')
folders_org = fp.readlines()

folders     = []
fp.close()

for i in range(len(folders_org)):
    if '#' in folders_org[i] or folders_org[i] == '\n':
        pass
    else:
        folders.append(folders_org[i].replace('\n',''))     

for path in folders:
    pwd = os.getcwd()
    os.chdir(path)
    print("process ", path)
    #s=System('POSCAR',fmt='poscar')
    ls=LabeledSystem('OUTCAR',fmt='outcar')
    ls.to_deepmd_raw('deepmd')
    ls.to_deepmd_npy('deepmd',set_size=1000)
    os.chdir(pwd)
    print("done ", path)

print("done")
