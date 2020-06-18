function mk_sne_cov,cosmo,sv

;Written Tom kitching, August 2008
;Adapted from Adam Amara's lensing cov code
;This routines calculates the error on m(z) for the supernovae calculation
;INPUT:  sv       (structure, version 0.11.1 and above)
;        cosmo    (structure, version 0.11.1 and above)
;OUTPUT: 

n_zbin=n_elements(sv.zbin_sn)  ; number of redshift bins

cov_z=dblarr(n_zbin)  ; covariance matrix

;error takes into accountdust, pec velocities, intrinsic scatter
for i=0,n_zbin-1 do begin
    cov_z(i)=((sv.sigmam^2.)+$
	((5*sv.sigmanu)/(2.303*(cosmo.const.c/1000.)*sv.zbin_sn[i]))^2.+$
	(sv.delm^2./sv.nbin_sn[i]))
endfor

; Store in structure
sne_cov={cov_z:cov_z}  ; binned quantities

return,sne_cov
end