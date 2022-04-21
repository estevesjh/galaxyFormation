#!/usr/bin/env python
import os
import pandas as pd
from astropy.table import Table
from astropy.io.fits import getdata
import astropy.io.ascii as at

class FileLocs(object):
    def __init__(self, dataset='sdss'):
        if dataset == 'sdss':
            self.data_loc = './data/catalogs/SDSS/'
            
            self.cluster = self.data_loc + 'groupCatalog_Yang_deCarvalho2017.csv'
            self.cluster_raw = self.data_loc + 'groupCatalog_Yang_deCarvalho2017.fits'
            self.cluster_sdss = self.data_loc + 'groupCatalog_Yang_deCarvalho2017_toSDSS.csv'
            
            self.galaxy_raw = self.data_loc + 'groups_deCarvalho2017_R200m_galaxies_final_flag_johnnyheq.csv'
            self.galaxy = self.data_loc + 'groups_deCarvalho2017_R200m_galaxies_final_flag_johnnyheq_pp.csv'
            self.galaxy_vl = self.data_loc + 'groups_deCarvalho2017_galaxies_final_flag_johnnyheq_pp_volumeLimited_v1.csv'
        
            self.morphology = self.data_loc + 'DL_morphology_SDSS_DS18.fit'

            self.datasets = {'cluster/raw': [self.cluster_raw,'fits'], 'cluster/main': [self.cluster,'csv'], 'cluster/sdss/': [self.cluster,'csv'],
                             'galaxy/raw': [self.galaxy_raw, 'csv'], 'galaxy/main': [self.galaxy_raw, 'csv'], 'galaxy/volumeLimited': [self.galaxy_raw, 'csv']}
            
        if dataset=='tng':
            self.data_loc = '../data/TNG/'
            self.hdf_fname_z0 = self.data_loc + 'TNG300-1_GalEvol_z0p00.hdf5'
            self.hdf_fname_p = self.data_loc + 'TNG300-1_GalEvol_z0p00_post.csv'
            self.cls_fname = self.data_loc + 'TNG300-1_GalEvol_z0p00_cluster_post.csv'
            
            #self.cat = load_files(self.hdf_fname_z0,key='Halos')
            # self.gal0 = at.read(self.hdf_fname_p)
            # self.cat = at.read(self.cls_fname)

    def load_catalogs(self, key):
        fname = self.datasets[key][0]
        dtype = self.datasets[key][1]
        
        if dtype == 'fits':
            data = Table(getdata(fname))
        if dtype == 'csv':
            data = at.read(fname)
        print(f'Loading Catalog: {fname}')
        return data
            
        
def load_files(fname,key='Halos'):
    pi = pd.read_hdf(fname, key=key)
    pi = Table.from_pandas(pi)
    return pi