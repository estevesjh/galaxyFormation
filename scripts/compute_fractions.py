import numpy as np
import pandas as pd
import os
from utils import gaussian_kde

class computeFraction:
    """This class provides a set of calculantions to compute the fraction and bootstrap errors
    """
    def __init__(self, name, path='../data/', key='all'):
        print('\nWelcome to our cluster enviromental effects tools')
        print('probablity: %s'%name)
        self.name = name
        self.path = path
        self.key  = key
        self.outfile_base = path+'outputs/%s/{name}_{var}.npy'%(key)
        check_dir(path+'outputs/%s/'%key)
    
    def add_probabilities(self, p, p_orbital, p_infall, p_interlopers, rn):
        self.data = {'p':p,'p_orbital':p_orbital,'p_infall':p_infall,'p_interlopers':p_interlopers,'rn':rn}
                        
    def run(self, xlabel, xvar, xbins, nBootStrap=10, write=False, is_field=False):
        # set bins
        nbins = xbins.size-1
        self.data[xlabel] = xvar
        
        # initiate dataframe
        df = pd.DataFrame(self.data)
        
        # initiate output arrays
        self.initiate_arrays(nsize=nbins, nBootStrap=nBootStrap)
        
        # compute fraction
        for i in range(nBootStrap):
            vec = df.sample(frac=1.0, replace=True).to_numpy().T

            xkeys, xmed = make_bins(vec[-1, :], xbins)
            
            frac1 = np.array([compute_fraction(vec[1, idx], vec[0, idx]) for idx in xkeys]).T
            frac2 = np.array([compute_fraction(vec[2, idx], vec[0, idx]) for idx in xkeys]).T
            frac3 = np.array([compute_fraction(vec[3, idx], vec[0, idx]) for idx in xkeys]).T
            frac4 = np.array([compute_fraction(vec[3, idx], vec[0, idx]) for idx in xkeys]).T
            frac12= np.array([compute_fraction(vec[1, idx]+vec[2, idx], vec[0, idx]) for idx in xkeys]).T
            
            frac4 = np.where(frac4<1.,np.median(frac4),frac4)
            
            if is_field:
                frac3 = frac4
                
            qfrac1 = quenching_fraction_excess(frac2, frac1)
            qfrac2 = quenching_fraction_excess(frac3, frac1)
            qfrac3 = quenching_fraction_excess(frac3, frac2)
            qfrac4 = quenching_fraction_excess(frac4, frac12)

            self.frac1[:, i] = frac1
            self.frac2[:, i] = frac2
            self.frac3[:, i] = frac3
            self.frac4[:, i] = frac4
        
            self.qfrac1[:, i] = qfrac1
            self.qfrac2[:, i] = qfrac2
            self.qfrac3[:, i] = qfrac3
            self.qfrac4[:, i] = qfrac4
        
        self.xmed = xmed
        
        if write:
            var = [self.xmed, self.frac1, self.frac2, self.frac3,\
                   self.qfrac1, self.qfrac2, self.qfrac3]
            #outfile = f'{self.path}/tmp/{xlabel}_{self.name}.npy'
            self.write(self.outfile_base, var, q1=16, q3=84)
            
        pass

    def run_kde(self, xlabel, xvar, xbins, bw=None, is_field=False, nBootStrap=10, write=False, radial_cut=False):
        # set bins
        nbins = xbins.size
        self.data[xlabel] = xvar
        
        # initiate dataframe
        df = pd.DataFrame(self.data)
        
        # initiate output arrays
        self.initiate_arrays(nsize=nbins, nBootStrap=nBootStrap)
        
        if radial_cut:
            inner_cut = xbins<=2.
            outer_cut = xbins>2.

        # compute fraction
        for i in range(nBootStrap):
            vec = df.sample(frac=1.0, replace=True).to_numpy().T
            frac1 = compute_fraction_kde(xbins, vec[-1], vec[1], vec[0], bw=bw)
            frac2 = compute_fraction_kde(xbins, vec[-1], vec[2], vec[0], bw=bw)
            frac3 = compute_fraction_kde(xbins, vec[-1], vec[3], vec[0], bw=bw)
            
            if radial_cut:
                frac12 = np.full((xbins.size,),np.nan)
                frac12[inner_cut] = frac1[inner_cut]
                frac4 = np.full((xbins.size,),np.nan)
                frac4[outer_cut] = frac3[outer_cut]
                
            else:
                vec1 = np.where(vec[-2]>2.,0.,vec[1])
                vec3 = np.where(vec[-2]<=2.,0.,vec[3])
                frac12 = compute_fraction_kde(xbins, vec[-1], vec1, vec[0], bw=bw)
                frac4 = compute_fraction_kde(xbins, vec[-1], vec3, vec[0], bw=bw)
            
            frac4_notnan= np.where(np.isnan(frac4),np.nanmedian(frac4),frac4)
            qfrac1 = quenching_fraction_excess(frac2, frac12)
            qfrac2 = quenching_fraction_excess(frac3, frac12)
            qfrac3 = quenching_fraction_excess(frac3, frac2)
            qfrac4 = quenching_fraction_excess(frac4_notnan, frac2)
                        
            self.frac1[:, i] = frac12
            self.frac2[:, i] = frac2
            self.frac3[:, i] = frac3
            self.frac4[:, i] = frac4
        
            self.qfrac1[:, i] = qfrac1
            self.qfrac2[:, i] = qfrac2
            self.qfrac3[:, i] = qfrac3
            self.qfrac4[:, i] = qfrac4
        
        self.xmed = xbins
        
        if write:
            var = [self.xmed, self.frac1, self.frac2, self.frac3, self.frac4,\
                   self.qfrac1, self.qfrac2, self.qfrac3, self.qfrac4]
            outfile = self.outfile_base.format(name=self.name,var=xlabel)
            self.write(outfile, var, q1=16, q3=84)
            
        pass

    def write(self, outfile, var, q1=16, q3=84):
        save_output(var, outfile, q1=q1, q3=q3)
    
    def initiate_arrays(self, nBootStrap=100, nsize=12):
        self.frac1 = np.full((nsize, nBootStrap), np.nan)
        self.frac2 = np.full((nsize, nBootStrap), np.nan)
        self.frac3 = np.full((nsize, nBootStrap), np.nan)
        self.frac4 = np.full((nsize, nBootStrap), np.nan)
        
        self.qfrac1 = np.full((nsize, nBootStrap), np.nan)
        self.qfrac2 = np.full((nsize, nBootStrap), np.nan)
        self.qfrac3 = np.full((nsize, nBootStrap), np.nan)
        self.qfrac4 = np.full((nsize, nBootStrap), np.nan)
        pass

def save_output(var,outfile,q1=16,q3=84):
    nfracs = len(var)-1
    nsize  = var[0].size
    out    = np.zeros((4*nfracs+1, nsize), dtype=float)

    out[0] = var[0]
    for i,x in enumerate(var[1:]):
        ii = 4*i+1
        out[ii] = np.nanmedian(x,axis=1)
        out[ii+1] = out[ii]-np.nanpercentile(x,q1,axis=1)
        out[ii+2] = np.nanpercentile(x,q3,axis=1)-out[ii]
    print(f'saved file: {outfile}')
    np.savetxt(outfile, out.T, fmt='%.5f')
    
def make_bins(x,xbins):
    indices = []
    xmd     = 0.5*(xbins[1:]+xbins[:-1])
    for xl,xh in zip(xbins[:-1],xbins[1:]):
        w, = np.where((x<=xh)&(x>xl))
        indices.append(w)
    return indices,xmd

def compute_fraction(prob1, prob2, eps=1e-3):
    N1, N2 = np.nansum(prob1*prob2), np.nansum(prob1)
    frac = N1/N2
    if N2<eps:
        frac = np.nan
    return frac

def compute_kde(x,weights,bw=0.1):
    Norm = np.nansum(weights)
    pdf = gaussian_kde(x, weights=weights, bw_method=bw)
    return pdf, Norm

def compute_fraction_kde(xvec,x,p1,p2,bw=None,eps=1e-3):
    if len(p1)>5:
        pdf1, N1 = compute_kde(x, p1, bw=bw)
        pdf2, N2 = compute_kde(x, p1*p2, bw=bw)
        denumerator = N1*pdf1(xvec)
        frac = N2*pdf2(xvec)/denumerator
        frac = np.where(denumerator<eps, np.nan, frac)
    else:
        frac = np.full((xvec.size,),np.nan)
    return frac


def quenching_fraction_excess(fq1,fq2):
    dfrac = fq2-fq1
    qfe = dfrac/(1-fq1)
    qfe = np.where(fq1==1,-1.,qfe)
    qfe = np.where(qfe>1, 1., qfe)
    qfe = np.where(qfe<-1, -1., qfe)
    return qfe

def check_dir(dir_name):
    if not os.path.isdir(dir_name):
        print(f'Creating directory: {dir_name}')
        os.makedirs(dir_name)
# q = computeFraction('quenching',path=fl.root+'outputs/',oufile_base=outfile_base)
# q.add_probabilities(qf,Po,Pi,Pn,Pf)
# q.run_kde('smass', mass, msbins, write=True, nBootStrap=100, bw=0.3)
# q.run_kde('radii', rn, rbins, write=True, nBootStrap=100, bw=0.2)