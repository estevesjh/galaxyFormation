import ast
import numpy as np

import sys
sys.path.append('./scripts')
from file_loc import FileLocs

from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology
from colossus.halo import concentration
from colossus.halo import mass_defs

import astropy.table as Table
from astropy import units as u
from astropy.constants import G, c

params = dict(H0 = 70, Om0 = 0.27, Ob0 = 0.0457, Tcmb0 = 2.7255, Neff = 3.046)
sigma8 = 0.82
ns = 0.96

astropy_cosmo = FlatLambdaCDM(**params)
astropy_cosmo.name = 'LCDM'
colossus_cosmo = cosmology.fromAstropy(astropy_cosmo, sigma8, ns, name = 'my_cosmo')

Msol = 1.98847e33
Mpc2cm = 3.086e+24
rad2deg= 180/np.pi
h=0.7
c_kms = c.value/1000

print(cosmology.getCurrent())
print('\n')

def AngularDistance(z):
    DA = float( (astropy_cosmo.luminosity_distance(z)/(1+z)**2)/u.Mpc ) # in Mpc
    return DA
AngularDistance = np.vectorize(AngularDistance)

#--- Critical universe density
def rhoc(z):
    try:
        rho_c = float(astropy_cosmo.critical_density(z)/(u.g/u.cm**3)) # in g/cm**3
    except:
        rho_c = [float(astropy_cosmo.critical_density(zi)/(u.g/u.cm**3)) for zi in z]
        rho_c = np.array(rho_c)
    
    rho_c = rho_c*(Mpc2cm**3)/Msol # in Msol/Mpc**3
    return rho_c

def convertM200toR200(M200,z, nc=200):
    ## M200 in solar masses
    ## R200 in Mpc
    rho = nc*rhoc(z)
    r200 = (3*M200/(4*np.pi/rho))**(1/3.)
    return r200/1e9

def convert_mass_defs(Mi,zi,mdef='200m'):
    ci = concentration.concentration(Mi, '200c', zi, model = 'diemer15')
    M200m, R200m, c200m = mass_defs.changeMassDefinition(Mi, ci, zi, '200c', mdef)
    return M200m, R200m, c200m

def vcirc(mass,redshift,mdef,cosmo):
    '''Calculate circular velocity in km/s for halos of mass M (Msun/h)'''
    rho_crit = cosmo.critical_density(redshift)
    if mdef[-1] == 'c':
        delta = int(mdef[:-1])
        rho = delta*rho_crit
    elif mdef[-1] == 'm':
        delta = int(mdef[:-1])
        rho = delta*rho_crit*cosmo.Om(redshift)
    else:
        raise RuntimeError("Not correct mass definition")
    v = np.sqrt(G*(np.pi*4*rho/3)**(1./3)*(mass*u.Msun)**(2./3))
    a = v.to(u.km/u.s)
    return a.value

if __name__ == "__main__":
    print('Starting Code')
    fl = FileLocs(dataset='sdss')
    data = fl.load_catalogs('cluster/raw')

    cid = np.array(data['Yang'])
    ra = np.array(data['_RAJ2000'])
    dec = np.array(data['_DEJ2000'])

    zcls = np.array(data['<z>'])
    r200c = np.array(data['R200'])/h
    m200c = np.array(10**data['logM200'])/h
    n200 = np.array(data['N200'])

    da     = AngularDistance(zcls)

    print('Converting to R200m')
    out = [convert_mass_defs(Mi,zi,'200m') for Mi,zi in zip(m200c,zcls)]
    m200m = np.array([line[0] for line in out])
    r200m = np.array([line[1]/1000. for line in out])
    c200m = np.array([line[2] for line in out])

    thetaR200c = (r200c*h/da)*(180/np.pi)
    thetaR200m = (r200m*h/da)*(180/np.pi)

    vcirc_c   = vcirc(m200c,zcls,'200c',astropy_cosmo)
    vcirc_m   = vcirc(m200m,zcls,'200m',astropy_cosmo)

    print('Converting to SDSS input')
    # Yang,RA,DEC,z,logM200,R200,thetaR200,N200
    columns = ['Yang','RA','DEC','redshift','logM200','thetaR200','thetaR200m']
    sdss = Table.QTable([cid,ra,dec,zcls,np.log10(m200c),thetaR200c,thetaR200m],names=columns)
    sdss.write(fl.cluster_sdss,format='csv',overwrite=True)
    print('SDSS table')
    sdss

    print('Creating main catalog table')
    columns = ['Yang','RA','DEC','redshift','N200','logM200c','logM200m','R200c','R200m','thetaR200','thetaR200m',"vcirc_c","vcirc_m"]
    cat = Table.QTable([cid,ra,dec,zcls,n200,np.log10(m200c),np.log10(m200m),r200c,r200m,thetaR200c,thetaR200m,vcirc_c,vcirc_m],names=columns)
    cat.write(fl.cluster,format='csv',overwrite=True)
    print('Cluster main table')
    cat
