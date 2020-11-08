import dpdata 
import numpy as np
import glob
import argparse
parser = argparse.ArgumentParser()
#l1=dpdata.System('/Users/jiedeng/Documents/ml/deepmd-kit/my_example/tail_3.dump')
parser.add_argument("--file","-f",type=str,default='dump.0',help="name of dump file,default dump.0")
parser.add_argument("--sigma","-s",type=float,help="sigma = kb*T")

args   = parser.parse_args()
l1=dpdata.System(args.file,fmt='lammps/dump')

l1.to_deepmd_npy('test')
l1.to_deepmd_raw('test')
sigma=args.sigma
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