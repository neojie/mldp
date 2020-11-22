#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 13:17:41 2020

@author: jiedeng
"""
import numpy as np
step=np.load('set.000/fparam.npy').shape[0]
print("{0} configs".format(step))
from subprocess import call
call("dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re7/pv_gpu.pb -d re7 -n {0}".format(step),shell=True)
call("dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re8/pv_gpu.pb -d re8 -n {0}".format(step),shell=True)
call("dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re5/pv_gpu.pb -d re5 -n {0}".format(step),shell=True)
call("dp test -m /u/home/j/jd848/project-lstixrud/pv_hf_copy/dp-train/extreme_filtered/re6/pv_gpu.pb -d re6 -n {0}".format(step),shell=True)
call("dp test -m /u/home/j/jd848/project-lstixrud/pv+hf/dp-train/model/scale2/extreme_filtered_re4/pv-gpu.pb -d re4 -n {0}".format(step),shell=True)
