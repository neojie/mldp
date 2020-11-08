import dpdata 
import numpy as np
import glob

#l1=dpdata.System('/Users/jiedeng/Documents/ml/deepmd-kit/my_example/tail_3.dump')
l1=dpdata.System('dump.0',fmt='lammps/dump')

l1.to_deepmd_npy('test')
l1.to_deepmd_raw('test')
sigma=0.344693
sets = glob.glob('test/*set*')
for seti in sets:
    print(seti)
    box=np.load(seti+'/box.npy')
    nsw = box.size//9
    print("the size of the file is: ", nsw )
    all_te = np.ones((nsw,))*sigma
    np.save( seti+'/fparam.npy', all_te)
print('done')


#ls=dpdata.System(lmp,fmt='lammps/dump')