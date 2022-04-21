import numpy as np

import sys
sys.path.append('./scripts')
from utils import check_non_valid_number
from file_loc import FileLocs
from compute_fractions import computeFraction

## set params
nBoot = 1000
run_name = 'volumeLimited'

## set bins
msbins = np.linspace(9.05,11.7,50)
rbins = np.arange(0.,3.0,0.1)+0.1

fl = FileLocs(dataset='sdss')
gal0 = fl.load_catalogs('galaxy/main')

mask  = (gal0['VLOS_MASK']).astype(bool)
mask &= gal0['VOLUME_LIM_MASK'].astype(bool)

gal = gal0[mask].copy()

print('Seting Variables')
print()
gid = np.array(gal['Yang'])

rn = np.array(gal['Rm'])
mass = np.array(gal['mass'])
# t_infall = np.array(gal['t_infall'])/1e9
morph_type = np.array(gal['TType'])
ssfr = np.array(gal['ssfr'])

# sfr classification
sf   = np.array(gal['SF']).astype(int)
qf   = (1-sf).astype(int)

# morphological classification
sp   = np.where(gal['TType'] > 0, 1, 0).astype(int)
ell  = np.where(gal['TType'] <=0, 1, 0).astype(int)
s0   = check_non_valid_number(gal['Pbulge'])
s0[np.isnan(s0)] = 0.

# b/t definition
bt = np.array(gal['BT'])
bt2 = np.where(bt>=0.5,1.,bt)
bt2 = np.where(bt<0.5,0.,bt2)

# dynamical probabilities
Pi   = np.array(gal['p_infall'])
Po   = np.array(gal['p_orbital'])
Pn   = np.array(gal['p_interlopers'])
m200 = np.array(gal['M200c'])

# mask
bt_mask = bt>=0.

q = computeFraction('quenching',path=fl.root,key=run_name)
q.add_probabilities(qf,Po,Pi,Pn,rn)
q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

q = computeFraction('elliptical',path=fl.root,key=run_name)
q.add_probabilities(ell,Po,Pi,Pn,rn)
q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

q = computeFraction('Bulge',path=fl.root,key=run_name)
q.add_probabilities(s0,Po,Pi,Pn,rn)
q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

# q = computeFraction('BT',path=fl.root,key=run_name)
# q.add_probabilities(bt[bt_mask],Po[bt_mask],Pi[bt_mask],Pn[bt_mask],rn[bt_mask])
# q.run_kde('smass', mass[bt_mask], msbins, write=True, nBootStrap=nBoot, bw=0.3)
# q.run_kde('radii', rn[bt_mask], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

q = computeFraction('BT_TH',path=fl.root,key=run_name)
q.add_probabilities(bt2[bt_mask],Po[bt_mask],Pi[bt_mask],Pn[bt_mask],rn[bt_mask])
q.run_kde('smass', mass[bt_mask], msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn[bt_mask], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)
