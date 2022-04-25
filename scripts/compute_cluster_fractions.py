# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at

import sys
sys.path.append('./scripts')
from utils import compute_fraction2, quenching_fraction_excess, chunks, check_non_valid_number
from file_loc import FileLocs

# select dataset
TNG = False
volumeLimited = True

dataset = 'sdss'
if TNG: dataset = 'tng'

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

labels_mpr = ['quenching','sf', 'bt']
mask_mpr = [qf, sf, bt2]
if not TNG:
    mask_mpr += [ell, sp, s0]
    labels_mpr += ['elliptical', 'spiral', 'bulge']

print('Compute Fractions')
print()
## Quenching Fraction
## i: infall, o: orbital, n: interlopers
keys = list(chunks(gid, cid))

# for each halo compute cluster, infall and interlopers fraction.
labels = labels_mpr#+labels_bpt
probs = mask_mpr#+mask_bpt

for l, p in zip(labels, probs):
    frac_o = np.array([compute_fraction2(Po[idx],p[idx]) for idx in keys]).T
    frac_i = np.array([compute_fraction2(Pi[idx],p[idx]) for idx in keys]).T
    frac_n = np.array([compute_fraction2(Pn[idx],p[idx]) for idx in keys]).T
    frac_f = np.array([compute_fraction2(Pf[idx],p[idx]) for idx in keys]).T
 
    qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
    qfrac_on = quenching_fraction_excess(frac_n, frac_o)
    qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    qfrac_if = quenching_fraction_excess(frac_f, frac_i)
    
    # save
    cat['fo_'+l] = frac_o[0]
    cat['fo_'+l+'_err'] = frac_o[1]
    cat['fi_'+l] = frac_i[0]
    cat['fi_'+l+'_err'] = frac_i[1]
    cat['fn_'+l] = frac_n[0]
    cat['fn_'+l+'_err'] = frac_n[1]
    cat['ff_'+l] = frac_f[0]
    cat['ff_'+l+'_err'] = frac_f[1]
    
    cat['qf1_'+l] = qfrac_oi[0]
    cat['qf1_'+l+'_err'] = qfrac_oi[1]
    cat['qf2_'+l] = qfrac_on[0]
    cat['qf2_'+l+'_err'] = qfrac_on[1]
    cat['qf3_'+l] = qfrac_in[0]
    cat['qf3_'+l+'_err'] = qfrac_in[1]
    cat['qf4_'+l] = qfrac_if[0]
    cat['qf4_'+l+'_err'] = qfrac_if[1]

print('Save file')
outname = fl.cluster_frac_th
if volumeLimited: outname = fl.cluster_frac_vl
cat.write(outname,overwrite=True)
print('Output name: %s'%outname)

print()
print('done!')