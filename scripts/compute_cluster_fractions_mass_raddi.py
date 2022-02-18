# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at
from utils import compute_fraction1, quenching_fraction_excess, chunks, check_non_valid_number

## stellar mass bin
msbins = np.arange(9.0,12.25,0.25)

## radii bins
rbins = np.arange(0.,3.25,0.25)

## free fall time bins
tbins = np.logspace(8.2,11.2,11)/1e9/2.   #infall time

## morph type
mrpbins = np.linspace(-3.,7.,11)

from file_loc import FileLocs
fl = FileLocs(dataset='tng')
# galaxy_file = fl.gal_fname1
# cluster_file = fl.cls_fname
outfile_base = fl.data_loc+'tmp/tng/{label}_{var}.npy'

print('--------Initial Files-------')
cat = fl.cat
gal0 = fl.gal0

print(gal0.colnames)
mask = np.abs(gal0['vlosn']) <= 3.
gal = gal0[mask].copy()

print('Seting Variables')
print()
cid = np.array(cat['HaloID'])
gid = np.array(gal['HostHaloID'])

rn = np.array(gal['Rn'])
mass = np.array(gal['mass'])
t_infall = np.array(gal['t_cross'])/1e9
#morph_type = np.array(gal['TType'])

# sfr classification
sf   = np.array(gal['SF']).astype(int)
qf   = (1-sf).astype(int)

# morphological classification
# sp   = np.where(gal['TType'] > 0, 1, 0).astype(int)
# ell  = np.where(gal['TType'] <=0, 1, 0).astype(int)
s0   = check_non_valid_number(gal['BT'])
# bulge= check_non_valid_number(gal['Pbulge'])
# disk = check_non_valid_number(gal['Pdisk'])
# bar  = check_non_valid_number(gal['PbarGZ2'])
# merger= check_non_valid_number(gal['Pmerg'])

# # bpt classification
# composite = (np.array(gal['bpt']) == 3).astype(int)
# lsf = ((np.array(gal['bpt']) == 1) | (np.array(gal['bpt']) == 2)).astype(int)
# liners = (np.array(gal['bpt']) == 5).astype(int)
# agn = (np.array(gal['bpt']) == 4).astype(int)
# unclass = (np.array(gal['bpt']) == -1).astype(int)

# orbital classification
Pi   = np.array(gal['pinfall'])
Po   = np.array(gal['porbital'])
Pn   = np.array(gal['pinterloper'])
m200 = np.array(gal['M200'])

# Po = np.where(gal['Rn']>=1.,0.,Po)

# quenched, spiral, SO, bulge, disco, barra, merger
# star forming, liners, AGN
# mass > 10. & mass z limited

labels_mpr = ['quenching','sf','bulge']
# labels_bpt = ['lsf', 'liners', 'agn', 'compos', 'unclas']

mask_mpr = [qf, sf, s0]

print('Compute Fractions')
print()
## Quenching Fraction
## i: infall, o: orbital, n: interlopers
from utils import compute_fraction, quenching_fraction_excess, save_output_matrix, make_bins

mskeys, mmed = make_bins(mass, msbins)
rkeys, rmed  = make_bins(rn, rbins)
tkeys, tmed = make_bins(t_infall, tbins)
# mrpkeys, mrpmed = make_bins(morph_type, mrpbins)

# for each halo compute cluster, infall and interlopers fraction.
labels = labels_mpr#+labels_bpt
probs = mask_mpr#+mask_bpt

print('Stellar Mass')
for l, p in zip(labels, probs):
    outfile = outfile_base.format(label=l, var='smass')
    print('Label: %s' % l)
    print('outfile: %s' % outfile)
    frac_o = np.array([compute_fraction1(Po[idx],p[idx]) for idx in mskeys]).T
    frac_i = np.array([compute_fraction1(Pi[idx],p[idx]) for idx in mskeys]).T
    frac_n = np.array([compute_fraction1(Pn[idx],p[idx]) for idx in mskeys]).T
 
    qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
    qfrac_on = quenching_fraction_excess(frac_n, frac_o)
    qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    
    # save
    save_output_matrix([mmed, frac_o, frac_i, frac_n, qfrac_oi, qfrac_on, qfrac_in],outfile)
print()

print('Radii')
for l, p in zip(labels, probs):
    outfile = outfile_base.format(label=l, var='radii')
    print('Label: %s' % l)
    print('outfile: %s' % outfile)
    frac_o = np.array([compute_fraction1(Po[idx],p[idx]) for idx in rkeys]).T
    frac_i = np.array([compute_fraction1(Pi[idx],p[idx]) for idx in rkeys]).T
    frac_n = np.array([compute_fraction1(Pn[idx],p[idx]) for idx in rkeys]).T
 
    qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
    qfrac_on = quenching_fraction_excess(frac_n, frac_o)
    qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    
    # save
    save_output_matrix([rmed, frac_o, frac_i, frac_n, qfrac_oi, qfrac_on, qfrac_in],outfile)
print()

print('Free Fall Time')
for l, p in zip(labels, probs):
    outfile = outfile_base.format(label=l, var='free_fall')
    print('Label: %s' % l)
    print('outfile: %s' % outfile)
    frac_o = np.array([compute_fraction1(Po[idx],p[idx]) for idx in tkeys]).T
    frac_i = np.array([compute_fraction1(Pi[idx],p[idx]) for idx in tkeys]).T
    frac_n = np.array([compute_fraction1(Pn[idx],p[idx]) for idx in tkeys]).T
 
    qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
    qfrac_on = quenching_fraction_excess(frac_n, frac_o)
    qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    
    # save
    save_output_matrix([tmed, frac_o, frac_i, frac_n, qfrac_oi, qfrac_on, qfrac_in], outfile)
print()

# print('Morphological Type')
# for l, p in zip(labels, probs):
#     outfile = outfile_base.format(label=l, var='morhp')
#     print('Label: %s' % l)
#     print('outfile: %s' % outfile)
#     frac_o = np.array([compute_fraction1(Po[idx],p[idx]) for idx in mrpkeys]).T
#     frac_i = np.array([compute_fraction1(Pi[idx],p[idx]) for idx in mrpkeys]).T
#     frac_n = np.array([compute_fraction1(Pn[idx],p[idx]) for idx in mrpkeys]).T
 
#     qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
#     qfrac_on = quenching_fraction_excess(frac_n, frac_o)
#     qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    
#     # save
#     save_output_matrix([mrpmed, frac_o, frac_i, frac_n, qfrac_oi, qfrac_on, qfrac_in], outfile)
# print()

# print('done!')