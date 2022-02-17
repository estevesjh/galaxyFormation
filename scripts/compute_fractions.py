import numpy as np
import pandas as pd

class computeFraction:
    """This class provides a set of calculantions to compute the fraction and bootstrap errors
    """
    def __init__(self, name, path='../data/'):
        print('Welcome to our cluster enviromental effects tools')
        print('probablity: %s'%name)
        self.name = name
        self.path = path
    
    def add_probabilities(self, p, p_orbital, p_infall, p_interlopers):
        self.data = {'p':p,'p_orbital':p_orbital,'p_infall':p_infall,'p_interlopers':p_interlopers}
                        
    def run(self, xlabel, xvar, xbins, nBootStrap=10, write=False):
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
            
            qfrac1 = quenching_fraction_excess(frac2, frac1)
            qfrac2 = quenching_fraction_excess(frac3, frac1)
            qfrac3 = quenching_fraction_excess(frac3, frac2)

            self.frac1[:, i] = frac1
            self.frac2[:, i] = frac2
            self.frac3[:, i] = frac3
        
            self.qfrac1[:, i] = qfrac1
            self.qfrac2[:, i] = qfrac2
            self.qfrac3[:, i] = qfrac3
        
        self.xmed = xmed
        
        if write:
            var = [self.xmed, self.frac1, self.frac2, self.frac3,\
                   self.qfrac1, self.qfrac2, self.qfrac3]
            outfile = f'{self.path}/tmp/{xlabel}_{self.name}.npy'
            self.write(outfile, var, q1=16, q3=84)
            
        pass
    
    def write(self, outfile, var, q1=16, q3=84):
        save_output(var, outfile, q1=q1, q3=q3)
    
    def initiate_arrays(self, nBootStrap=100, nsize=12):
        self.frac1 = np.full((nsize, nBootStrap), np.nan)
        self.frac2 = np.full((nsize, nBootStrap), np.nan)
        self.frac3 = np.full((nsize, nBootStrap), np.nan)
        
        self.qfrac1 = np.full((nsize, nBootStrap), np.nan)
        self.qfrac2 = np.full((nsize, nBootStrap), np.nan)
        self.qfrac3 = np.full((nsize, nBootStrap), np.nan)
        pass

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

def quenching_fraction_excess(fq1,fq2):
    dfrac = fq2-fq1
    qfe = dfrac/(1-fq1)
    qfe = np.where(fq1==1,np.nan,qfe)
    qfe = np.where(qfe>1, 1., qfe)
    qfe = np.where(qfe<-1, -1., qfe)
    return qfe

def save_output(var,outfile,q1=16,q3=84):
    nfracs = len(var)-1
    nsize  = var[0].size
    out    = np.zeros((3*nfracs+1, nsize), dtype=float)

    out[0] = var[0]
    for i,x in enumerate(var[1:]):
        ii = 3*i+1
        out[ii] = np.nanmedian(x,axis=1)
        out[ii+1] = out[ii]-np.nanpercentile(x,q1,axis=1)
        out[ii+2] = np.nanpercentile(x,q3,axis=1)-out[ii]
    np.savetxt(outfile, out.T, fmt='%.5f')