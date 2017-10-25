import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

path = "MBAR_estimates_all_samples.csv"

df = pd.read_csv(path, sep=';')

eps = np.array(df['epsilon values'])
rmin_half = np.array(df['rmin_half values'])

eps_ = np.array(sorted(list(set(eps))))
rmin_half_ = np.array(sorted(list(set(rmin_half))))

empty_matrix = lambda : np.zeros((len(eps_), len(rmin_half_)))

C_p = empty_matrix()
vol = empty_matrix()

dC_p = empty_matrix()
dvol = empty_matrix()

n_eff = empty_matrix()

eps_dict = {}
for i in range(len(eps_)):
    eps_dict[eps_[i]] = i
rmin_half_dict = {}
for i in range(len(rmin_half_)):
    rmin_half_dict[rmin_half_[i]] = i

for ind in range(len(eps)):
    i = eps_dict[eps[ind]]
    j = rmin_half_dict[rmin_half[ind]]
    
    C_p[i,j] = df['C_p_expect (J/mol/K)'][ind]
    vol[i,j] = df['Vol_expect (mL/mol)'][ind]
    
    dC_p[i,j] = df['dC_p_expect (J/mol/K)'][ind]
    dvol[i,j] = df['dVol_expect (mL/mol)'][ind]
    
    n_eff[i,j] = df['N_eff'][ind]

def add_axis_labels_and_colorbar():
    plt.xlabel(r"$r_m / 2$")
    plt.ylabel(r"$\epsilon$")
    plt.colorbar()

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.title('C_p_expect (J/mol/K)')
plt.contourf(rmin_half_, eps_, C_p, 30, cmap = 'viridis')
plt.plot(1.908,0.1094,'x',color='red')
add_axis_labels_and_colorbar()

plt.subplot(1,2,2)
plt.title('dC_p_expect (J/mol/K)')
plt.contourf(rmin_half_, eps_, dC_p, 30, cmap = 'viridis')
plt.plot(1.908,0.1094,'x',color='orange')
add_axis_labels_and_colorbar()

plt.tight_layout()
plt.savefig("c_p.png", dpi=300)

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.title('Vol_expect (mL/mol)')
plt.contourf(rmin_half_, eps_, vol, 30, cmap = 'viridis')
plt.plot(1.908,0.1094,'x',color='yellow')
add_axis_labels_and_colorbar()

plt.subplot(1,2,2)
plt.title('dVol_expect (mL/mol)')
plt.contourf(rmin_half_, eps_, dvol, 30, cmap = 'viridis')
plt.plot(1.908,0.1094,'x',color='orange')
add_axis_labels_and_colorbar()
plt.tight_layout()
plt.savefig("vol.png", dpi=300)

plt.title('N_eff')
plt.contourf(rmin_half_, eps_, n_eff, 30, cmap = 'viridis')
plt.plot(1.908,0.1094,'x',color='purple')
add_axis_labels_and_colorbar()
plt.savefig("n_eff.png", dpi=300)

plt.plot(n_eff,dvol,'bo')
plt.xlabel('# of effective samples')
plt.ylabel('d(Volume) (mL/mol)')
plt.savefig('dvol_vs_neff.png', dpi=300)

plt.plot(n_eff,dC_p,'bo')
plt.xlabel('# of effective samples')
plt.ylabel('d(C_p) (J/mol/K)')
plt.savefig('dCp_vs_neff.png', dpi=300)
