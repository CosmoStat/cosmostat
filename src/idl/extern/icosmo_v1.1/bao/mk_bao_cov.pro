;pro mk_bao_cov,bao_rad, bao_trans, cosmo,survey=survey, deltay, deltayp
function mk_bao_cov,cosmo,bao,survey
  
; Jul 08 - modified by AA to be consisted with v0.11
;Written by Anais Rassat, July 2008
;Calculates the transverse and radial errors for the BAO scale.
;NEED MORE INFO
;INPUT:
;OUTPUT: x = Delta y, xp = Delta yp
  
;OPTIONAL INPUT:
; bao: type of BAO Calculation
;     0: Spectroscopic Tangential
;     1: Spectroscopic Radial
;     2: Photometric Tangential

;TO DO:
; At the moment doesn't have input for photometric calculation


  
;Determine survey bins
  deltayp=mk_cov_blake(cosmo,survey,1) ; Spectroscopic Radial
  deltay=mk_cov_blake(cosmo,survey,0)  ; Spectroscopic Tangential
  
;Blake et al 2006 outputs x = delta y / y *100, and for the fisher
;matrix I want deltay
  deltay = deltay/100.d0*bao.y_bin
  deltayp = deltayp/100.d0*bao.yprime_bin
  
  bao_cov={delta_trans: deltay, delta_rad:deltayp}
return,bao_cov
end
