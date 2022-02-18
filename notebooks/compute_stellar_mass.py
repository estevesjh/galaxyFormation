# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at
from astropy.table import Table

from utils import compute_fraction2, quenching_fraction_excess, chunks, check_non_valid_number

from file_loc import FileLocs
fl = FileLocs()
galaxy_file = fl.gal_fname1
cluster_file = fl.cls_fname

cluster_file_3 = cluster_file.split('.csv')[0]+'_smass.csv'

print('--------Initial Files-------')
print('Cluster File : %s' % cluster_file)
print('Galaxy File : %s' % galaxy_file)
print('Outfile: %s' % cluster_file_3)
print()

# load catalogs
cat0 = at.read(cluster_file)
gal0 = at.read(galaxy_file)

mask = np.abs(gal0['vlosn']) <= 3.
gal = gal0[mask].copy()

print('Seting Variables')
print()
cid = np.array(cat0['Yang'])
gid = np.array(gal['Yang'])

mass = np.array(gal0['mass'])
mass_err = np.array((gal0['mass_p84']-gal0['mass_p16'])/2.)

sfr = np.array(gal0['sfr'])
sfr_err = np.array((gal0['sfr_p84']-gal0['sfr_p16'])/2.)

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

# Po = np.where(gal['Rn']>=1,0.,Po)
# Pi = np.where(gal['Rn']<=1,0.,Pi)
Pn = np.where(gal['Rn']<=1,0.,Pn)

# quenched, spiral, SO, bulge, disco, barra, merger
# star forming, liners, AGN
# mass > 10. & mass z limited

labels = ['quenching', 'sf', 'elliptical', 'spiral']
mask = [qf, sf, ell, sp]

print('Compute Fractions')
print()
## Quenching Fraction
## i: infall, o: orbital, n: interlopers
cat = Table(cat0[['Yang','z','N200','logM200']])
keys = list(chunks(gid, cid))

def compute_weighted_sum(weigths,xlog,xlog_err):
    x = 10**xlog
    x_err = xlog_err*x
    x_sum = np.nansum(weigths*x)
    weights_n = weigths/np.nansum(weigths)
    x_sum_err = np.sqrt(np.nansum((x_err*weights_n)**2))
    x_sum_err = x_sum_err/x_sum
    return np.array([np.log10(x_sum),x_sum_err])

# for each halo compute cluster, infall and stellar mass content
smass_o = np.array([compute_weighted_sum(Po[idx],mass[idx],mass_err[idx]) for idx in keys]).T
smass_i = np.array([compute_weighted_sum(Pi[idx],mass[idx],mass_err[idx]) for idx in keys]).T
sfr_o = np.array([compute_weighted_sum(Po[idx],sfr[idx],sfr_err[idx]) for idx in keys]).T
sfr_i = np.array([compute_weighted_sum(Pi[idx],sfr[idx],sfr_err[idx]) for idx in keys]).T

cat['smasso'] = smass_o[0]
cat['smasso_err'] = smass_o[1]
cat['smassi'] = smass_i[0]
cat['smassi_err'] = smass_i[1]

cat['sfro'] = sfr_o[0]
cat['sfro_err'] = sfr_o[1]
cat['sfri'] = sfr_i[0]
cat['sfri_err'] = sfr_i[1]

for l, p in zip(labels, mask):
    smass_o = np.array([compute_weighted_sum(p[idx]*Po[idx],mass[idx],mass_err[idx]) for idx in keys]).T
    smass_i = np.array([compute_weighted_sum(p[idx]*Pi[idx],mass[idx],mass_err[idx]) for idx in keys]).T
    sfr_o = np.array([compute_weighted_sum(p[idx]*Po[idx],sfr[idx],sfr_err[idx]) for idx in keys]).T
    sfr_i = np.array([compute_weighted_sum(p[idx]*Pi[idx],sfr[idx],sfr_err[idx]) for idx in keys]).T
 
    # save
    cat['smasso_'+l] = smass_o[0]
    cat['smasso_'+l+'_err'] = smass_o[1]
    cat['smassi_'+l] = smass_i[0]
    cat['smassi_'+l+'_err'] = smass_i[1]
    
    cat['sfro_'+l] = sfr_o[0]
    cat['sfro_'+l+'_err'] = sfr_o[1]
    cat['sfri_'+l] = sfr_i[0]
    cat['sfri_'+l+'_err'] = sfr_i[1]

print('Save file')        
cat.write(cluster_file_3,overwrite=True)

print()
print('done!')