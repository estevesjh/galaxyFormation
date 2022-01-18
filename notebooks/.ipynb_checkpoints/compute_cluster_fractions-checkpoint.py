# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at

from utils import compute_fraction2, quenching_fraction_excess, chunks, check_non_valid_number

from file_loc import FileLocs
fl = FileLocs()
galaxy_file = fl.gal_fname1
cluster_file = fl.cls_fname

cluster_file_2 = cluster_file.split('.csv')[0]+'_frac.csv'

print('--------Initial Files-------')
print('Cluster File : %s' % cluster_file)
print('Galaxy File : %s' % galaxy_file)
print('Outfile: %s' % cluster_file_2)
print()

# load catalogs
cat = at.read(cluster_file)
gal0 = at.read(galaxy_file)

mask = np.abs(gal0['vlosn']) <= 3.
gal = gal0[mask].copy()

print('Seting Variables')
print()
cid = np.array(cat['Yang'])
gid = np.array(gal['Yang'])

# sfr classification
sf   = np.array(gal['SF']).astype(int)
qf   = (1-sf).astype(int)

# morphological classification
sp   = np.where(gal['TType'] > 0, 1, 0).astype(int)
ell  = np.where(gal['TType'] <=0, 1, 0).astype(int)
s0   = check_non_valid_number(gal['PS0'])
bulge= check_non_valid_number(gal['Pbulge'])
disk = check_non_valid_number(gal['Pdisk'])
bar  = check_non_valid_number(gal['PbarGZ2'])
merger= check_non_valid_number(gal['Pmerg'])

# bpt classification
composite = (np.array(gal['bpt']) == 3).astype(int)
lsf = ((np.array(gal['bpt']) == 1) | (np.array(gal['bpt']) == 2)).astype(int)
liners = (np.array(gal['bpt']) == 5).astype(int)
agn = (np.array(gal['bpt']) == 4).astype(int)
unclass = (np.array(gal['bpt']) == -1).astype(int)

# orbital classification
Pi   = np.array(gal['pinfall'])
Po   = np.array(gal['porbital'])
Pn   = np.array(gal['pinterloper'])
m200 = np.array(gal['M200'])

Po = np.where(gal['Rn']>=1,0.,Po)
Pi = np.where(gal['Rn']<=1,0.,Pi)

# quenched, spiral, SO, bulge, disco, barra, merger
# star forming, liners, AGN
# mass > 10. & mass z limited

labels_mpr = ['quenching', 'sf', 'elliptical', 'spiral', 's0', 'bulge', 'disk', 'bar', 'merger']
labels_bpt = ['lsf', 'liners', 'agn', 'compos', 'unclas']

mask_mpr = [qf, sf, ell, sp, s0, bulge, disk, bar, merger]
mask_bpt = [lsf, liners, agn, composite, unclass]

print('Compute Fractions')
print()
## Quenching Fraction
## i: infall, o: orbital, n: interlopers

keys = list(chunks(gid, cid))

# for each halo compute cluster, infall and interlopers fraction.
labels = labels_mpr+labels_bpt
probs = mask_mpr+mask_bpt

for l, p in zip(labels, probs):
    frac_o = np.array([compute_fraction2(Po[idx],p[idx]) for idx in keys]).T
    frac_i = np.array([compute_fraction2(Pi[idx],p[idx]) for idx in keys]).T
    frac_n = np.array([compute_fraction2(Pn[idx],p[idx]) for idx in keys]).T
 
    qfrac_oi = quenching_fraction_excess(frac_i, frac_o)
    qfrac_on = quenching_fraction_excess(frac_n, frac_o)
    qfrac_in = quenching_fraction_excess(frac_n, frac_i)
    
    # save
    cat['fo_'+l] = frac_o[0]
    cat['fo_'+l+'_err'] = frac_o[1]
    cat['fi_'+l] = frac_i[0]
    cat['fi_'+l+'_err'] = frac_i[1]
    cat['fn_'+l] = frac_n[0]
    cat['fn_'+l+'_err'] = frac_n[1]
    
    cat['qf1_'+l] = qfrac_oi[0]
    cat['qf1_'+l+'_err'] = qfrac_oi[1]
    cat['qf2_'+l] = qfrac_on[0]
    cat['qf2_'+l+'_err'] = qfrac_on[1]
    cat['qf3_'+l] = qfrac_in[0]
    cat['qf3_'+l+'_err'] = qfrac_in[1]

print('Save file')        
cat.write(cluster_file_2,overwrite=True)

print()
print('done!')