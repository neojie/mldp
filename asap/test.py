#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 18:45:26 2021

@author: jiedeng
"""

test = ASAPXYZ('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap_test_for_size/ASAP-desc.xyz')

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e',use_atomic_desc=True)[0].shape

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e')

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e',use_atomic_desc=True,species_name = 8)




test1 = ASAPXYZ('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap_test_for_size/test1.xyz')

test1.get_descriptors('SOAP-n6-l6-c6.0-g0.44',use_atomic_desc=True)

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e')

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e',use_atomic_desc=True,species_name = 8)


test2 = ASAPXYZ('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap_test_for_size/test2.xyz')

test2.get_descriptors('SOAP-n6-l6-c6.0-g0.44',use_atomic_desc=True)

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e')

test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e',use_atomic_desc=True,species_name = 8)


peratom1 = ASAPXYZ('/Users/jiedeng/Documents/tmp/jd848/project_folder/pv+hf/3k/solid1/r3-3k/recal/asap_test_for_size/per_atom.xyz')

peratom1.get_descriptors('SOAP-n6-l6-c6.0-g0.44',use_atomic_desc=True) # somehow 'SOAP-n6-l6-c6.0-g0.44-e' doe snot work


t10_no_e = test1.get_descriptors('SOAP-n6-l6-c6.0-g0.44',use_atomic_desc=True)[0][0] # first frame
t0_with_e = test.get_descriptors('SOAP-n6-l6-c6.0-g0.44-e',use_atomic_desc=True)[0][0] # first frame

p10_mg  = peratom1.get_atomic_descriptors('SOAP-n6-l6-c6.0-g0.44',species_name = 12)[:32]
p10_si  = peratom1.get_atomic_descriptors('SOAP-n6-l6-c6.0-g0.44',species_name = 14)[:32]
p10_o  = peratom1.get_atomic_descriptors('SOAP-n6-l6-c6.0-g0.44',species_name = 8)[:96]





### no -e param
tmp =np.concatenate((p10_mg,p10_si,p10_o))
plt.plot(tmp.mean(axis=0))
plt.plot(t10_no_e)

### with -e param


tmp =np.concatenate((p10_o.mean(axis=0), p10_mg.mean(axis=0),p10_si.mean(axis=0)))
plt.plot(tmp)
plt.plot(t0_with_e)
