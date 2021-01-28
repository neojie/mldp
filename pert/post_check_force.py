#import dpdata as dp
#from dpdata import LabeledSystem
#outcar = '/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/8k/r3-3100/recal/OUTCAR'
#ls = LabeledSystem(outcar,fmt='outcar')
#tmp = ls[0]

import numpy as np

force=np.load('set.000/force.npy')
force1=np.load('set.001/force.npy')
print('set0')
print(force.max(axis=1))
print(force.min(axis=1))
print(min(force.min(axis=1)), max(force.max(axis=1)))
print(np.argmin(force.min(axis=1)), np.argmax(force.max(axis=1)))
print('set1')
print(force1.max(axis=1))
print(force1.min(axis=1))
print(min(force1.min(axis=1)), max(force1.max(axis=1)))
print(np.argmin(force1.min(axis=1)), np.argmax(force1.max(axis=1)))


#force=np.load('/Users/jiedeng/GD/papers/pv3_crystallization/post_nn/generate_new_poscar/pairs/untitled folder/sisi/deepmd/set.000/force.npy')

#max_f = []
#for f in force:
#    max_f.append(np.max(np.max(f)))
#print(max_f)    
#plt.plot(max_f)
tmp  = '/Users/jiedeng/GD/papers/pv3_crystallization/post_nn/exsolution_pert/u.project.ESS.lstixrud.jd848.pv+hf.dp-train.lmp_run.6k.rp5.160-cpu.pert.10k_good_p3.recal/deepmd_all'

tmp2 = '/Users/jiedeng/GD/papers/pv3_crystallization/post_nn/exsolution_pert/u.home.j.jd848.project-lstixrud.metad.3rd.recal/deepmd/'
ls = LabeledSystem(tmp2,fmt='deepmd/npy')

from dpdata import LabeledSystem
import numpy as np
ls = LabeledSystem('deepmd_ttr2',fmt='deepmd/npy')
idx = list(range(len(ls)))
idx.remove(14)
ls2 = ls.sub_system(idx)
ls2.to_deepmd_npy('test')
fparam = np.load('deepmd_ttr2/set.001/fparam.npy')
np.save('test/set.000/fparam.npy',fparam[idx])

cp -r deepmd_ttr2 deepmd_ttr2_tmp
rm -r deepmd_ttr2_tmp/set.000

from dpdata import LabeledSystem
import numpy as np
ls = LabeledSystem('deepmd_ttr2_tmp',fmt='deepmd/npy')
idx = list(range(len(ls)))
idx.remove(10)
ls2 = ls.sub_system(idx)
ls2.to_deepmd_npy('test')
fparam = np.load('deepmd_ttr2_tmp/set.000/fparam.npy')
np.save('test/set.000/fparam.npy',fparam[idx])
