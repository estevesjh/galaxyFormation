#!/usr/bin/env python
import os
import pandas as pd
from astropy.table import Table
import astropy.io.ascii as at

class FileLocs(object):
    def __init__(self, dataset='sdss'):
        if dataset=='sdss':
            self.data_loc = '../data/'
            self.cls_fname   = self.data_loc + 'groupCatalog_Yang_deCarvalho2017.csv'
            self.gal_fname0  = self.data_loc + 'groups_deCarvalho2017_galaxies_final_flag_johnnyheq.csv'
            self.gal_fname1  = self.data_loc + 'groups_deCarvalho2017_galaxies_final_flag_johnnyheq_volumeLimited_v1.csv'
            self.morp_fname = self.data_loc + 'DL_morphology_SDSS_DS18.fit'

            # load catalogs
            self.cat = at.read(self.cls_fname)
            self.gal0 = at.read(self.gal_fname1)

        if dataset=='tng':
            self.data_loc = '../data/'
            self.hdf_fname_z0 = self.data_loc + 'TNG300-1_GalEvol_z0p00.hdf5'
            self.hdf_fname_p = self.data_loc + 'TNG300-1_GalEvol_z0p00_post.csv'
            self.cls_fname = self.data_loc + 'TNG300-1_GalEvol_z0p00_cluster_post.csv'
            
            #self.cat = load_files(self.hdf_fname_z0,key='Halos')
            self.gal0 = at.read(self.hdf_fname_p)
            self.cat = at.read(self.cls_fname)

def load_files(fname,key='Halos'):
    pi = pd.read_hdf(fname, key=key)
    pi = Table.from_pandas(pi)
    return pi