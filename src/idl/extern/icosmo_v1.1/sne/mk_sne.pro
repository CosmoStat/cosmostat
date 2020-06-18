function mk_sne,fid,cosmo,sv

; Sep 08 - Modified by AnR: check_sv added
;Merger of mk_sne made from mk_bao, written by Tom Kitching Aug 2008
;Calculates binned dl and converts to magnitudes

;Find out if survey is also defined for supernovae calculation:
check_sv, sv, 'sne'


if not keyword_set(sv) then stop,'Error: suvey structure has not been included'

zmid_bin=sv.zbin_sn; mid point of the bins 

;convert to H0 independant magnitudes
factor = 1000.d/cosmo.const.c

dl_bin = interpol(cosmo.evol.dl, cosmo.evol.z, zmid_bin)
z_bin  = interpol(cosmo.evol.z, cosmo.evol.z, zmid_bin)

mz_bin =5.*alog10(factor*dl_bin)

mz =5.*alog10(factor*cosmo.evol.dl)

M0 =5.*alog10(fid.cosmo.h*100)

sne={mz_bin:mz_bin,z_bin:z_bin,mz:mz,M0:M0}

  return, sne
end

  
