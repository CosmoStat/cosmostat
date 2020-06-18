function mk_pk_grid,fid,cosmo,plotit=plotit

; Jul 08 - Modified by AA - added k_ran to structure
; July 07 - Written by A. Refregier
;
; PURPOSE: Compute the nonlinear and linear normalised matter power
; spectrum P(k,z) on a grid of wavenumber k and redshifts z. This
; grid will is to be appended to the cosmo strucure. Different
; prescriptions and fitting functions are supported for the linear power
; spectrum and for non-linear corrections.
; INPUT: fid: fiducial parameter structure from set_fiducial.pro
;        cosmo: intermediate cosmological parameter structure from
;               mk_cosmo.pro
; OPTIONAL INPUT: plotit: plot P(k,z) for debugging
; OUTPUT: pk: structure containing the linear and non-linea power spectrum
;         cosmo.cost.norma: gets updated with correct normalisation
; TO DO: add examples and add fiducial parameter used

; declarations
n_z=fid.calc.nz_crs  ; number of z -values to consider 
z=xgen(fid.calc.ran_z(0),fid.calc.ran_z(1),n=n_z,/double)
n_k=fid.calc.n_k
pk_l_arr=dblarr(n_k,n_z)        ; linear power spectrum P_l(k,z) array
pk_nl_arr=dblarr(n_k,n_z)       ; nonlinear power spectrum P_nl(k,z) array

; set up old calculation
;del2_old=mk_del2_new(fid,cosmo,0.)
;n_k_old=n_elements(del2_old.k)
;pk_l_arr_old=dblarr(n_k_old,n_z)           ; linear power spectrum P_l(k,z) array
;pk_nl_arr_old=dblarr(n_k_old,n_z)          ; nonlinear power spectrum P_nl(k,z) array

; compute power spectrum at each z
;stop
for i=0,n_z-1 do begin
  ; compute linear and non-linear power spectrum and normalise 
  ; the former in the process
  pk_nl_zi=mk_pk_nl(fid,cosmo,z=z(i))  

  ; fill in arrays
  pk_l_arr(*,i)=pk_nl_zi.pk_l
  pk_nl_arr(*,i)=pk_nl_zi.pk

  ; old calculation
;  del2_old=mk_del2_new(fid,cosmo,z(i))
;  pk_l_arr_old(*,i)=del2_old.pk_l
;  pk_nl_arr_old(*,i)=del2_old.pk

endfor 

; reconstructing the range used to make k vector
; In the current implimentation (Jul 2008) this is 
; exactly the same k_ran in fiducial structure from set_fiducial
k_pk=pk_nl_zi.k
del_logk=alog(k_pk(1:*))-alog(k_pk)
del_logk=del_logk(uniq(float(del_logk)))
if (n_elements(del_logk) ne 1) then stop ;!!!del is changing!!!
k_rantemp=[k_pk(0),exp(alog(k_pk(0))+del_logk*n_elements(k_pk))] 
if (abs(total((fid.calc.k_ran-k_rantemp)/fid.calc.k_ran)) gt 1.e-8) then stop ;!!! at the moment k_ran is only set by set_fiducial so the k_ran for the pk calculation and fiducial structure should match!!!
k_ran=fid.calc.k_ran


; store result in a structure
pk={k:pk_nl_zi.k,z:z,pk:pk_nl_arr,pk_l:pk_l_arr,k_ran:k_ran}
;pk_old={k:del2_old.k,pk:pk_nl_arr_old,pk_l:pk_l_arr_old}

; plot power spectra if requested
if keyword_set(plotit) then begin
  plot,[0],[0],xran=[1e-4,1e6],yran=[1e-15,1e5],/xlog,/ylog,/nodata,$   
    xtitle='k [h Mpc^-1]',ytitle='P(k)'
  tek_color
  for i=0,n_z-1 do begin
    oplot,pk.k,pk.pk_l(*,i)
    oplot,pk.k,pk.pk(*,i)
;    oplot,pk_old.k,pk_old.pk_l(*,i),col=2
;    oplot,pk_old.k,pk_old.pk(*,i),col=2
  endfor
endif

return,pk
end



