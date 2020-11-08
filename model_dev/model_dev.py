import numpy as np
import glob
es_dirs=glob.glob('*.e.out')
fs_dirs=glob.glob('*.f.out')
vs_dirs=glob.glob('*.v.out')
es = []
fs = []
vs = []
for dire in es_dirs:
    tmp = np.loadtxt(dire)
    es.append(tmp[:,1])
es = np.array(es)

for dire in fs_dirs:
    tmp = np.loadtxt(dire)
    tmp = tmp[:,3:6]
    yy =tmp.reshape((100,160,3))
    zz = yy.reshape((100,480))
    fs.append(zz)    
fs = np.array(fs)
for dire in vs_dirs:
    tmp = np.loadtxt(dire)
    tmp = tmp[:,9:]
    vs.append(tmp)   
vs = np.array(vs)


var_e = np.var(es,axis=0)
std_e = np.std(es,axis=0)
var_f = np.var(fs,axis=0);max_var_f = np.max(var_f,axis=1)
std_f = np.std(fs,axis=0);max_std_f = np.max(std_f,axis=1)
var_v = np.var(vs,axis=0);max_var_v = np.max(var_v,axis=1)
std_v = np.std(vs,axis=0);max_std_v = np.max(std_v,axis=1)

import matplotlib.pyplot as plt
plt.plot(std_e,label='e')
plt.plot(max_std_f,label='f')
plt.plot(max_std_v,label='v')
plt.legend() 
plt.show()
plt.savefig('test.png')




"""
dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re7/pv_gpu.pb -d re7
dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re8/pv_gpu.pb -d re8
dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re5/pv_gpu.pb -d re5
dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re6/pv_gpu.pb -d re6
dp test -m /u/home/j/jd848/project-lstixrud/pv+hf/dp-train/model/scale2/extreme_filtered_re4/pv-gpu.pb -d re4
"""