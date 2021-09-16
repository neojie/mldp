import matplotlib.pyplot as plt
import pandas as pd
import os
out_pd = pd.read_csv('out/out_pd')
target='out'
fig, ax = plt.subplots(1,3,figsize=(9,3))
ax[1].plot(out_pd['T'],out_pd['P(GPa)'],'o')
ax[0].plot(out_pd['T'],out_pd['vol'],'o')
ax[2].plot(out_pd['T'],out_pd['E(eV)'],'o')
ax[0].set_xlabel('T'); ax[0].set_ylabel('P (GPa)')
ax[1].set_xlabel('T'); ax[0].set_ylabel('Vol (A3)')
ax[2].set_xlabel('T'); ax[0].set_ylabel('E (eV)')
fig.savefig(os.path.join(target,'pev.png'),bbox_inches='tight')    

fig, ax = plt.subplots(1,3,figsize=(9,3))
ax[0].scatter(out_pd['sigma(eV)'],out_pd['Est_eq'],c=out_pd['vol'])
ax[1].scatter(out_pd['sigma(eV)'],out_pd['Fst_eq'],c=out_pd['vol'])
f3=ax[2].scatter(out_pd['sigma(eV)'],out_pd['Vst_eq'],c=out_pd['vol'])
plt.colorbar(f3,ax[2])
ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Est (GPa)')
ax[1].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Fst (A3)')
ax[2].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Vst (eV)')
fig.savefig(os.path.join(target,'simga_vs_std.png'),bbox_inches='tight')  

fig, ax = plt.subplots(1,3,figsize=(9,3))
ax[0].scatter(out_pd['sigma(eV)'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
ax[1].scatter(out_pd['sigma(eV)'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
f3=ax[2].scatter(out_pd['sigma(eV)'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
plt.colorbar(f3,ax[2])
ax[0].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Eaerr/std (GPa)')
ax[1].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Fst/std (A3)')
ax[2].set_xlabel('sigma(eV)'); ax[0].set_ylabel('Vst/std (eV)')
fig.savefig(os.path.join(target,'simga_vs_err/std.png'),bbox_inches='tight')  


fig, ax = plt.subplots(1,3,figsize=(9,3))
ax[0].scatter(out_pd['Ndata'],out_pd['Eaerr_tr']/(out_pd['Est_eq']/out_pd['N']),c=out_pd['vol'])
ax[1].scatter(out_pd['Ndata'],out_pd['Ferr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
f3=ax[2].scatter(out_pd['Ndata'],out_pd['Verr_tr']/(out_pd['Fst_eq']/out_pd['N']),c=out_pd['vol'])
plt.colorbar(f3,ax[2])
ax[0].set_xlabel('Ndata'); ax[0].set_ylabel('Eaerr/std (GPa)')
ax[1].set_xlabel('Ndata'); ax[0].set_ylabel('Fst/std (A3)')
ax[2].set_xlabel('Ndata'); ax[0].set_ylabel('Vst/std (eV)')
fig.savefig(os.path.join(target,'ndata_vs_err/std.png'),bbox_inches='tight') 

