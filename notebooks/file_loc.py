#!/usr/bin/env python
import os

class FileLocs(object):
    def __init__(self, machine='fnal'):
        self.data_loc = '../data/'
        self.cls_fname   = self.data_loc + 'groupCatalog_Yang_deCarvalho2017.csv'
        self.gal_fname0  = self.data_loc + 'groupCatalog_Yang_deCarvalho2017_galaxy.csv'