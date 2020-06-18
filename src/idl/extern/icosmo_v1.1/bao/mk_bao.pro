function mk_bao,cosmo,survey
;Sep 08 - Modified by AnR: check_sv added
;Merger of mk_bao_rad and mk_bao_trans, written by Anais Rassat July
;2008
; ***************// HELP BEGIN //**************
; PURPOSE: Calculates the BAO Radial and Transverse scales, both
;          binned and unbinned.
;          For the radial mode: calculates y = r(z)/s 
;          For the transverse mode: calculates y' = c/H(z)/sh 
;          Note: s is sound horizon
; INPUT: cosmo: structure output by cosmo= mk_cosmo(fid)
;        survey: structure output by structure=mk_survey(fid)
; OPTIONAL INPUT: None
; OUTPUT: bao: structure containing:
;             yprime: transverse BAO scale as a function of redshift
;             y:      radial BAO scale as a function of redshift
;             z:      redshift range
;             yprime_bin: binned transverse BAO scale as a function of redshift
;             y_bin:      binned radial BAO scale as a function of redshift
;             z_bin:      binned redshift range
; OPTIONAL OUTPUT: None
; NOTE: Output can be plotted using plt_bao
;---
; Example 
;            > fid = set_fiducial()
;            > cosmo = mk_cosmo(fid)
;            > sv = mk_survey(fid, 'sv1')
;            > bao = mk_bao(cosmo, sv)
; ***************// HELP END //**************

;Check if survey is designed for BAO calculation: 
check_sv, survey, 'bao'

if not keyword_set(survey) then stop,'Error: suvey structure has not been included'
;if (strlowcase(survey.dndztype) ne 'hist') then stop,'Only valid for histogram inputs of pz.'

zmid_bin=(survey.z_min+survey.z_max)/2.d0; mid point of the bins 
                                ;!! NOT SURE THIS IS WHAT WE WANT IF
                                ;WE ARE NOT IN HIST MODE FOR P(Z) !!

yprime = 1.d0/cosmo.evol.dzdr/cosmo.const.sh ; Radial: y' = c/H(z)/sh 
yprime_bin = interpol(yprime, cosmo.evol.z, zmid_bin)

y  = cosmo.evol.sk*cosmo.const.r0/cosmo.const.sh ; Tangential: y = r(z)/s
y_bin = interpol(y, cosmo.evol.z, zmid_bin)

;create structure
bao={yprime:yprime,y:y,z:cosmo.evol.z,yprime_bin:yprime_bin,y_bin:y_bin,z_bin:zmid_bin}

  return, bao
end

  
