----------------------------------------------

Data from SDSS + Yang et al. 2005 + deCarvalho et al. 2017

----------------------------------------------

HOW DATA WAS CONSTRUCTED:

	1. Cluster data (from groupCatalog_Yang_deCarvalho2017.fits) 
    is uploaded from the deCarvalho et al. 2017 paper. 

	2. Using the query (getClusterGalaxies.sql) on the CasJobs (SDSS) platform 
    I uploaded all the galaxies within 2xR200 in the photometric (photoObjAll) 
    and spectroscopic (specObjAll).

    3. The derived spectral quantites are taken from two tables, galSpecExtra (MPA-JHU SDSS DR7)
    and eBOSS/SDSS firefly (DR14/DR16).


HOW TO READ THIS DATA:

    The csv and fits can be read with astropy

    Under construstion

---------------------------------------

Spectral Catalog properties

---------------------------------------
---------------------------------------
A: MPA-JHU SDSS DR7 

The catalog have their estimations of stellar mass, star formation rate (SFR), specific SFR and BPT class.
Also, the 16 and 84 percentiles of each quantity. 

For more information: https://www.sdss.org/dr16/spectro/galaxy_mpajhu/

The files and tables above also include a number of galaxy parameters derived by the MPA-JHU group available.
•	BPT classification: We supply emission line classifications based on the BPT diagram. 
    Galaxies are divided into “Star Forming”, “Composite”, “AGN”, “Low S/N Star Forming”, 
    “Low S/N AGN”, and “Unclassifiable” categories.
•	Stellar Mass: Stellar masses are calculated using the Bayesian methodology and model grids described 
    in Kauffmann et al. (2003). The spectra are measured through a 3 arcsec aperture, and therefore 
    do not represent the entire galaxy. We therefore base our model on the ugriz galaxy photometry 
    alone (rather than the spectral indices Dn(4000) and Hδ used by Kauffmann et al. 2003). 
    We have corrected the photometry for the small contribution due to nebular emission using the spectra. 
    We estimate the stellar mass within the SDSS spectroscopic fiber aperture using fiber magnitudes and 
    the total stellar mass using model magnitudes. A Kroupa (2001) initial mass function is assumed. 
    We output the stellar mass corresponding to the median and 2.5%, 16%, 84%, 97.5% of the probability 
    distribution function.
•	Star Formation Rate: Star formation rates (SFRs) are computed within the galaxy fiber aperture 
    using the nebular emission lines as described in Brinchmann et al. (2004). SFRs outside of the fiber 
    are estimated using the galaxy photometry following Salim et al. (2007). For AGN and galaxies with 
    weak emission lines, SFRs are estimated from the photometry. We report both the fiber SFR and 
    the total SFR at the median and 2.5%, 16%, 84%, 97.5% of the probability distribution function.
•	Specific SFR: The Specific SFR (SFR divided by the stellar mass) has been calculated by combining 
    the SFR and stellar mass likelihood distributions as outlined in Appendix A of Brinchmann et al. (2004). 
    We report both the fiber and the total specific SFR at the median and 2.5%, 16%, 84%, 97.5% of 
    the probability distribution function.

---------------------------------------
---------------------------------------
B: eBOSS SDSS firefly

The firefly code is spectral energy density (SED) fitting code. So their estimated quantities are based on
the spectrum information not on the photometry information as in A. We use their mass and age estimations.
I could not find their star formation history estimation altough they describe in their paper (Comparat et al. 2017)

For more information: https://www.sdss.org/dr16/spectro/eboss-firefly-value-added-catalog/#DR16Products



