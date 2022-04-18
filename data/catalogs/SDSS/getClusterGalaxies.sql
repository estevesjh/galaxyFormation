SELECT
  m.Yang, n.objid, s.specobjid, n.distance, o.ra AS ra, o.dec AS dec, o.petroMag_r as mr_petro,
  o.dered_u as mu, o.dered_g as mg, o.dered_r as mr, o.dered_i as mi, o.dered_z as mz,
  o.modelMagErr_u as mu_Err, o.modelMagErr_g as mg_err, o.modelMagErr_r as mr_err, o.modelMagErr_i as mi_err, o.modelMagErr_z as mz_err,
   c.z as photo_z, c.zErr as photo_zErr, 
  ISNULL(s.z,-99) as z,ISNULL(s.zErr,-99) as zErr, o.fracDeV_i, 
  mpa.bptclass as bpt, ISNULL(mpa.lgm_tot_p50,-99) as mass, ISNULL(mpa.lgm_tot_p16,-99) as mass_p16, ISNULL(mpa.lgm_tot_p84,-99) as mass_p84,
  ISNULL(mpa.sfr_tot_p50,-99) as sfr, ISNULL(mpa.sfr_tot_p16,-99) as sfr_p16, ISNULL(mpa.sfr_tot_p84,-99) as sfr_p84, 
  ISNULL(mpa.specsfr_tot_p50,-99) as ssfr, ISNULL(mpa.specsfr_tot_p16,-99) as ssfr_p16, ISNULL(mpa.specsfr_tot_p84,-99) as ssfr_p84, 
  ISNULL(fire.Chabrier_MILES_total_mass,-99) as mass_fire, ISNULL(fire.Chabrier_MILES_total_mass_low_1sig,-99) as mass_fire_p16, ISNULL(fire.Chabrier_MILES_total_mass_up_1sig,-99) as mass_fire_p84,
  ISNULL(fire.Chabrier_MILES_age_lightW,-99) as age_fire, ISNULL(fire.Chabrier_MILES_age_lightW_low_1sig,-99) as age_fire_p16, ISNULL(fire.Chabrier_MILES_age_lightW_up_1sig,-99) as age_fire_p84,
  o.score as score
   into mydb.groups_deCarvalho2017_R200m_galaxies from MyDB.deCarvalho2017_R200m AS m
     CROSS APPLY dbo.fGetNearbyObjAllEq( m.RA, m.DEC, 2.1*60*m.thetaR200m) AS n
      JOIN Galaxy AS o ON n.objid = o.objid
      left outer join SpecObj as s on o.objID = s.bestObjID
      join Photoz as c on c.objID = o.objID
      left outer join galSpecExtra as mpa on s.specObjID = mpa.specObjID
      left outer join sdssEbossFirefly as fire on s.specObjID = fire.SPECOBJID
WHERE
  o.petroMag_r < 17.80 AND
  ( (o.dered_g - o.dered_r) BETWEEN -1 AND 4) AND
  (dbo.fPhotoFlags('SATUR_CENTER')) != 0 AND 
  (dbo.fPhotoFlags('BRIGHT')) != 0 AND 
  (dbo.fPhotoFlags('NODEBLEND')) != 0 AND 
  (dbo.fPhotoFlags('DEBLEND_TOO_MANY_PEAKS')) != 0
  AND ((o.flags_r & 0x10000000) != 0)
  
  AND ((o.flags_r & 0x8100000c00a0) = 0)