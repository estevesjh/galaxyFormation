#!/usr/bin/env python
import os

class FileLocs(object):
    def __init__(self, machine='fnal'):
        self.data_loc = '../data/'
        self.cls_fname   = self.data_loc + 'groupCatalog_Yang_deCarvalho2017.csv'
        self.gal_fname0  = self.data_loc + 'groups_deCarvalho2017_galaxies_final_flag_johnnyheq.csv'
        #self.gal_fname0  = self.data_loc + 'groupCatalog_Yang_deCarvalho2017_galaxy.csv'
        #self.gal_fname1  = self.data_loc + 'groupCatalog_Yang_deCarvalho2017_galaxy_v1.csv'
        self.gal_fname1  = self.data_loc + 'groups_deCarvalho2017_galaxies_final_flag_johnnyheq_volumeLimited_v1.csv'
        self.morp_fname = self.data_loc + 'DL_morphology_SDSS_DS18.fit'