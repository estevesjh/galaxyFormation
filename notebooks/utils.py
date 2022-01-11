# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at

## stellar mass bin
msbins = np.arange(9.0,12.25,0.25)

## radii bins
rbins = np.arange(0.,3.25,0.25)

def compute_fraction(prob,eps=1e-6):
    N1, N2 = np.nansum(prob), prob[~np.isnan(prob)].size
    frac = N1/N2
    frac_err = frac*np.sqrt(1/(N1+eps)+1/(N2+eps))
    return np.array([frac,frac_err])

def quenching_fraction_excess(fq1,fq2):
    dfrac = fq2[0]-fq1[0]
    qfe = dfrac/(1-fq1[0])
    qfe_err = qfe*np.sqrt((np.sqrt(fq2[1]**2+fq1[1]**2)/dfrac)**2 + (fq1[1]/(1-fq1[0]))**2)
    return np.array([qfe,qfe_err])

def make_bins(x,xbins):
    indices = []
    xmd     = 0.5*(xbins[1:]+xbins[:-1])
    for xl,xh in zip(xbins[:-1],xbins[1:]):
        w, = np.where((x<=xh)&(x>xl))
        indices.append(w)
    return indices,xmd

def chunks(ids1, ids2):
    """Yield successive n-sized chunks from data"""
    for id in ids2:
        w, = np.where( ids1==id )
        yield w
        
def save_output_matrix(var,outfile):
    nout = len(var)
    nsize = var[0].size
    out = np.zeros((2*nout, nsize) ,dtype=float)
    
    out[0] = var[0]
    for i,x in enumerate(var[1:]):
        ii = 2*i+1
        out[ii] = x[0]
        out[ii+1] = x[1]
    np.savetxt(outfile, out.T, fmt='%.5f')
    
def check_non_valid_number(x):
    w,  = np.where(x < 0)
    x[w] = np.nan
    return np.array(x)

label1 = ['Cluster', 'Infall', 'Interlopers']
label2 = ['Cluster+Infall', 'Cluster+Interlopers', 'Infall+Interlopers']

def plot_fraction_pannel(label, type):
    x = np.loadtx(file_base.format(label, type))
    
    fig = plt.figure(figsize=(12,4))

    plt.subplot(1, 2, 1)
    for i in range(3):
        ii = 2*i + 1
        plt.errorbar(x[0], x[ii], yerr=x[ii+1], label=label1[i], fmt='o')
    
    plt.legend(fontsize=12,loc=2)
    plt.xlabel(r'Log($M_\star/M_{\odot}$)',fontsize=16)
    plt.title('Quenching Fraction',fontsize=16)
    plt.ylim(-0.02,1.05)

    plt.subplot(1, 2, 2)
    for i in range(3, 6):
        ii = 2*i + 1
        plt.errorbar(x[0], x[ii], yerr=x[ii+1], label=label2[i], fmt='o')
    plt.legend(fontsize=12,loc=2)
    plt.xlabel(r'Log($M_\star/M_{\odot}$)',fontsize=18)
    plt.title('Quenching Fraction Excess',fontsize=16)

    plt.ylim(-0.02,1.05)