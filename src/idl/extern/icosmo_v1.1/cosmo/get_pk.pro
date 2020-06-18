function get_pk,cosmo,z=z,plotit=plotit

; July 08 - Written by A. Refregier
; July 08 - Modified by A. Rassat, removed fiducial as an input.

; PURPOSE: Extract the matter power spectrum by interpolating
; the cosmo.pk grid P(k,z) onto a given redshift z. This is done for
; both the linear and nonlinear power spectra.
; INPUT:  cosmo: intermediate cosmological parameter structure from
;               mk_cosmo.pro
; OPTIONAL INPUT: z: redshift (default: z=0)
;                 plotit: plot P(k,z) for debugging
; OUTPUT: pk: structure containing the linear and non-linea power spectrum

; define location arrays
if not keyword_set(z) then z=0.
n_z=n_elements(cosmo.pk.z)   ; number of z vaues on the grid
if z lt cosmo.pk.z(0) or z gt cosmo.pk.z(n_z-1) then $
  stop,'ERROR: get_pk.pro: z outside grid range'

k=cosmo.pk.k        ; wavenumber k on the grd [h Mpc^-1]
n_k=n_elements(k)   ; number of k values on the grid
ik=indgen(n_k)      ; k indices: output k values are 
                    ; identical as that in the grid     
iz=interpol(findgen(n_z),cosmo.pk.z,z) ; index value corresponding to requested
                                      ; z on the grid
iz=replicate(iz,n_k)

; interpolate grid
pk_l=interpolate(cosmo.pk.pk_l,ik,iz)  ; linear power spectrum
pk_nl=interpolate(cosmo.pk.pk,ik,iz)  ; linear power spectrum

; compute Delta^2
del2_l=k^3*pk_l/(2.d*!dpi^2.d)
del2_nl=k^3*pk_nl/(2.d*!dpi^2.d)

; direct/old computation
;del2_old=mk_del2_new(fid,cosmo,z)
;pk_direct=mk_pk_nl(fid,cosmo,mk_pk_lin(fid,cosmo,z=z),z)

; store result in a structure
pk={z:z,k:k,pk:pk_nl,del2:del2_nl,pk_l:pk_l,del2_l:del2_l}

; plot power spectra
if keyword_set(plotit) then begin
  plot,[0],[0],xran=[1e-4,1e6],yran=[1e-15,1e5],/xlog,/ylog,/nodata,$   
    xtitle='k [h Mpc^-1]',ytitle='P(k)'
;  tek_color
  oplot,cosmo.pk.k,pk_l
  oplot,cosmo.pk.k,pk_nl
;  oplot,del2_old.k,del2_old.pk_l,color=2
;  oplot,del2_old.k,del2_old.pk,color=2
;  oplot,pk_direct.k,pk_direct.pk_l,color=3
;  oplot,pk_direct.k,pk_direct.pk,color=3
endif


return,pk
end



