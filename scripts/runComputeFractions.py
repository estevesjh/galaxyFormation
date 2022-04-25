from re import M
import numpy as np

#python scripts/runComputeFractions.py > log.out &
import sys
sys.path.append('./scripts')
from utils import check_non_valid_number
from file_loc import FileLocs
from compute_fractions import computeFraction

## set params
##############################################
TNG = False
volumeLimited = True

nBoot = 1000
##############################################
msbins = np.linspace(10.05,11.7,50)
rbins = np.arange(0.,3.0,0.1)+0.1

dataset = 'sdss'
run_name = 'magTH'

if TNG: 
    dataset = 'tng'
    
if volumeLimited: 
    run_name = 'volumeLimited'
    msbins = np.linspace(9.05,11.7,50)

## set bins
# msbins = np.linspace(10.35,11.7,50)
# msbins = np.linspace(9.05,11.7,50)

fl = FileLocs(dataset=dataset)
cat = fl.load_catalogs('cluster/main')
gal0 = fl.load_catalogs('galaxy/main')

mask  = (gal0['VLOS_MASK']).astype(bool)
mask &= gal0['mass']>=9.

if volumeLimited:
    mask &= gal0['VOLUME_LIM_MASK'].astype(bool)
else:
    mask &= gal0['MAG_TH_MASK'].astype(bool)

gal = gal0[mask].copy()

print('Seting Variables')
print()
print('Seting Variables')
print()
cid = np.array(cat['Yang'])
gid = np.array(gal['Yang'])

rn = np.array(gal['Rm'])
mass = np.array(gal['mass'])
ssfr = np.array(gal['ssfr'])

# sfr classification
sf   = np.array(gal['SF']).astype(int)
qf   = (1-sf).astype(int)

# morphological classification
if not TNG:
    morph_type = np.array(gal['TType'])
    sp   = np.where(gal['TType'] > 0, 1, 0).astype(int)
    ell  = np.where(gal['TType'] <=0, 1, 0).astype(int)
    s0   = check_non_valid_number(gal['Pbulge'])
    s0[np.isnan(s0)] = 0.
    
    ps0 = np.array(gal['PS0'])
    ps0 = np.where(ps0<0.,0.,ps0)
    ell_pure = np.where(ell,1-ps0,0.)

# b/t definition
bt = np.array(gal['BT'])
bt2 = np.where(bt>=0.5,1.,bt)
bt2 = np.where(bt<0.5,0.,bt2)

# dynamical probabilities
Pi   = np.array(gal['p_infall'])
Po   = np.array(gal['p_orbital'])
Pn   = np.array(gal['p_interlopers'])
m200 = np.array(gal['M200c'])

Pf   = np.where(rn<2.,0.,Pn)
Po   = np.where(rn>2.,0.,Po)

# mask
bt_mask = bt>=0.

fl.root = fl.root+dataset.upper()+'/'
q = computeFraction('quenching',path=fl.root,key=run_name)
q.add_probabilities(qf,Po,Pi,Pn,rn)
q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

q = computeFraction('BT',path=fl.root,key=run_name)
q.add_probabilities(bt2,Po,Pi,Pn,rn)
q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

if not TNG:
    q = computeFraction('elliptical',path=fl.root,key=run_name)
    q.add_probabilities(ell,Po,Pi,Pn,rn)
    q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
    q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

    q = computeFraction('elliptical_pure',path=fl.root,key=run_name)
    q.add_probabilities(ell_pure,Po,Pi,Pn,rn)
    q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
    q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

    q = computeFraction('S0',path=fl.root,key=run_name)
    q.add_probabilities(ps0,Po,Pi,Pn,rn) 
    q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
    q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

    q = computeFraction('Bulge',path=fl.root,key=run_name)
    q.add_probabilities(s0,Po,Pi,Pn,rn) 
    q.run_kde('smass', mass, msbins, write=True, nBootStrap=nBoot, bw=0.3)
    q.run_kde('radii', rn, rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)
    
    pass

## mass thresholds
if volumeLimited:
    mCuts = np.array([9.,10.05, 10.53, 13.])
    mass_keys = [np.where((mass>=mlow)&(mass<mhig))[0] for mlow,mhig in zip(mCuts[:-1],mCuts[1:])]
    mass_label = ['mLow','mMean','mHigh']

    for idx, run_name in zip(mass_keys,mass_label):
        q = computeFraction('quenching',path=fl.root,key=run_name)
        q.add_probabilities(qf[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
        q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

        q = computeFraction('BT',path=fl.root,key=run_name)
        q.add_probabilities(bt2[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
        q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

    if not TNG:
        for idx, run_name in zip(mass_keys,mass_label):
            q = computeFraction('elliptical',path=fl.root,key=run_name)
            q.add_probabilities(ell[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
            q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

            q = computeFraction('elliptical_pure',path=fl.root,key=run_name)
            q.add_probabilities(ell_pure[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
            q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)

            q = computeFraction('S0',path=fl.root,key=run_name)
            q.add_probabilities(ell_pure[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
            q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)
            
            q = computeFraction('Bulge',path=fl.root,key=run_name)
            q.add_probabilities(ell_pure[idx],Po[idx],Pi[idx],Pn[idx],rn[idx])
            q.run_kde('radii', rn[idx], rbins, write=True, nBootStrap=nBoot, bw=0.2,radial_cut=True)
