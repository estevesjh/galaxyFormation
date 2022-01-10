#############################################################################
                          
		Morphological catalogue for SDSS galaxies
                     
			current version: v1.0, Feb 2018

#############################################################################

The morphologies have been obtained by applying Deep Learning models to SDSS-DR7 RGB cutouts, as explained in Dominguez Sánchez+18.

The catalogue contains 670,722 rows, one per galaxy.
The meaning of the columns are  (see also Table 3 from DS+18):

=====================================
Name                 Meaning
=====================================

 1  - dr7objid     =  SDSS-DR7 ID
 2  - galcount     =  Meert+15 ID
 3  - P_disk       =  Probability of showing disk/features (vs. being smooth)
 4  - P_edge_on    =  Prob. of being edge-on 
 5  - P_bar_GZ2    =  Prob. of having bar signature (trained with GZ2 catalogue)
 6  - P_bar_Nair10 =  Prob. of having bar signature (trained with Nair+10 catalogue)
 7  - P_merg       =  Prob. of being merger/projected pairs
 8  - P_bulge      =  Prob. of having a dominant/obvious bulge (vs. no bulge)
 9  - P_cigar      =  Prob. of having cigar shape (vs. round shape)
 10 - TType        =  T-Type
 11 - P_S0         =  Prob of being S0 (vs. pure Ell)


Notes:
———————————————————————————
7  - P_merg;  proxy to clustered galaxies or projected pairs. Need for visual inspection to select true on-going mergers

10 - TType values range = [-3.3, 8.0]

	* ETG: TType <= 0
	* LTG: Type > 0

11 - P_S0 separates between Ell and S0. Only meaningful for T-Type <= 0




