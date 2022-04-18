----------------------------------------------

Data for Johnny Galaxy Evolution study

----------------------------------------------


HOW TO READ THIS DATA:

	The HDF5 files were written using the pandas python package.
	You can quickly load them and turn them into 2D numpy arrays
	using the following commands:

	*****

	import pandas as pd
	import numpy as np

	Halos    = pd.read_hdf('PATH TO HDF5 FILE', key = 'Halos').to_numpy()
	Subhalos = pd.read_hdf('PATH TO HDF5 FILE', key = 'Subhalos').to_numpy()


	*****

	Each of 'Halos' and 'Subhalos' will be a 2D numpy array 
	that you can index into using slicing. The column index of each
	quantity is given below in the square brackets '[i]'

	If you leave them as pandas dataframes you can also quickly access
	an individual column and/or turn it into a numpy array via:

	******

	Pandas_Series = Halos[quantity_name]
	Numpy_array   = Halos[quantity_name].values

	
	******
	


NOTES:
	1. TNG runs under a standard Planck15 cosmology. Using "Planck15"
	   in Astropy, or Colossus, provides the required parameters.

	2. Some fields may contain the value np.NaN.
	   This is a placeholder value/label where a quantity was not computed
	   for various reasons.

---------------------------------------

Halo Catalog properties

	Contains Halos with M200c > 10^14 Msun
---------------------------------------

HaloID[0]:
	Unique ID for the halo at THIS redshift. 
	The HaloID is NOT unique across all snapshots.

M200c[1]:
	Total mass (in log10(Msun/h) units) enclosed in 
	sphere within which density is 200 times critical density

R200c[2]:
	Radius (in physical kpc) at which enclosed mass 
	is 200 times critical density

---------------------------------------

Subhalo catalog properties

	Only contains subhalos with MStargal > 10^9 Msun,
	and subhalos within 3*R200c (projected radius) of each halo.
	Include interlopers as well. Galaxy quantities are only computed 
	for subhalos with abs(vlos) < 4*V_c, where V_c is the 
	circular velocity of host halo. R_proj and V_los are computed
	for all subhalos.
---------------------------------------

SubhaloID[0]:
	ID of this subhalo in this snapshot.

MStargal[1]: 
	SUBFIND total stellar mass of subhalo. In log10(Msun/h).
	Is the total mass of all gravitationally bound particles
	to this subhalo.

HostHaloID[2]:
	ID of this subhalo's host halo. A single SubhaloID
	can appear multiple times with different HostHaloID
	each time given the projection we do.

R_Proj[3]:
	The projected radius of the subhalo from the host halo
	(or the BCG, since central subhalo and host halo share
	same position/velocity). Is given in units of R200c so
	it is dimensionless. Central subhalos always have R_proj = 0.

v_los[4]:
	The line of sight velocities (always along x-axis) of galaxies,
	with Hubble flow velocities added. In units km/s

V200c_host[5]:
	The circular velocity of the host halo, V200c = sqrt(G*M200c/R200c). In units km/s

t_infall[6]: 
	The time (in Gyrs) when the subhalo cross R200c (in 3D space)
	of the host halo. If a subhalo never crossed R200c (i.e. it is
	always inside or outside R200c) it is assigned a np.NaN infall time.

Mstar_bulge[7]: 
	The sum of stellar mass associated with the spheroidal bulge. (In log10(Msun) units)
	Identified as particles with j_z/j < 0.7, where j is the total angular
	momentum of particle, and j_z is the component aligned with the total angular
	momentum vector of all star particles. We also apply a correction factor
	of M_true = M_measure/0.85. See arxiv:1904.12860 for details. Only particles
	within r < 3*R_1/2 are considered, whre R_1/2 the stellar half-mass radius.

S_over_T[8]:
	The ratio of bulge mass to total mass, obtained as Mstar_bulge/Mstar_tot.
	
	NOTE: Mstar_tot is not the same as Mstargal. The former is computed within 
	R < 3*R_1/2 whereas the latter uses all particles bound to subhalo.	

SFR[9]: 
	The total instantaneous star formation rate of all gas cells in ths subhalo.
	Note that individual gas cells with SFR < 1e-3 are assigned SFR = 0 due to
	resolution issues. In units Msun/yr

Stellar_age[10]:
	The stellar mass-weighted age of the stars in a subhalo. In units of Gyrs.