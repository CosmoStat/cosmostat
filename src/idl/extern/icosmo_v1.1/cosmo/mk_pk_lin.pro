function mk_pk_lin,fid,cosmo,z=z

; Aug 08 - modified by AA - make sure sigma_8 calc has sufficient
;          k_ran and k resolution
; Jul 08 - created by A. Refregier from mk_del2_lin.pro
;
; PURPOSE: compute the linear matter power spectrum P(k) at a given
; redshift z.This is done using either the BBKS (as described in
; Peacock and Dodds, 1996) or the Eisenstein & Hu (1998) transfer functions
; INPUT: fid: fiducial parameter structure from set_fid.pro
;        cosmo: cosmological parameter structure from mk_cosmo.pro
; OPTIONAL INPUT: z: redshift (default: z=0)
; OUTPUT: pk_lin: power spectrum structure containing:
;            k: comoving wave number in h-corrected physical units
;                i.e. k=k_comoving/R_0/h   [h Mpc^-1]
;            pk_l: linear power spectrum [h^-3 Mpc^3]
;            del2_l: variance contribution per log k defined as
;                    Delta^2(k)=k^3*P(k)/(2pi^2) [1]

;stop

; declarations:
n_k=fid.calc.n_k             ;number of k values to consider
k=xgen(fid.calc.k_ran(0),fid.calc.k_ran(1),n=n_k,/log,/double)
                                  ; k values to consider [h Mpc^-1]
dlnk=(alog(fid.calc.k_ran(1))-alog(fid.calc.k_ran(0)))/n_k

k_ran_sigma8=[.001d,100.0d] ; minimum range of k-values to compute sigma8 [h Mpc^-1]
dlnk_sigma8=0.01d ;minimum log k sampling to compute sigma8
if not keyword_set(z) then z=0.

; compute linear power spectrum at z=0, unnormalised
pk_l=k^cosmo.const.n*$
       tk(k,cosmo.const,fit_tk=fid.calc.fit_tk)^2.d*(2.*!pi^2) 
pk_lin={z:z,k:k,pk_l:pk_l}     ; temporary structure

; normalise if necessary
if cosmo.const.norma eq -1.d then begin
   if cosmo.const.sigma8 eq -1. then begin ; CMB normalisation
      stop,'CMB normalisation not yet implemented'
   endif else begin             ; sigma8 normalisation
      
      if fid.calc.k_ran(0) lt k_ran_sigma8(0) and $
         fid.calc.k_ran(1) gt k_ran_sigma8(1) and $
         dlnk lt dlnk_sigma8 then pk_lin_temp=pk_lin $
      else begin
         k_ran_temp=[0.d,0.d]
         k_ran_temp(0)=min([fid.calc.k_ran(0),k_ran_sigma8(0)])
         k_ran_temp(1)=max([fid.calc.k_ran(1),k_ran_sigma8(1)])
         dlnk_temp=min([dlnk,dlnk_sigma8])
         n_k_temp=ceil((alog(k_ran_temp(1)) - alog(k_ran_temp(0)))/dlnk_temp)
         if (n_k gt 1000000L) then stop,'Error, something is wrong with the range and resolution in k that you have chosen'
         k_temp=xgen(k_ran_temp(0),k_ran_temp(1),n=n_k_temp,/log,/double)
         pk_l_temp=k_temp^cosmo.const.n*$
                   tk(k_temp,cosmo.const,fit_tk=fid.calc.fit_tk)^2.d*$
                   (2.d *!dpi^2)          
         pk_lin_temp={z:z,k:k_temp,pk_l:pk_l_temp} ; temporary structure 
      endelse
                                ; compute sigma8
;      stop
      sig8=sigma8(pk_lin_temp)  
   endelse
   cosmo.const.norma=(cosmo.const.sigma8/sig8)^2
endif
pk_l=pk_l*cosmo.const.norma

; multiply by growth factor if z<>0
if z ne 0. then pk_l=pk_l*interpol(cosmo.evol.d,cosmo.evol.z,z)^2.d

; compute Delta^2(k), the contribution to the variance per log k
del2_l=k^3*pk_l/(2.d*!dpi^2.d)

; construct structure
pk_lin={z:z,k:k,pk_l:pk_l,del2_l:del2_l}

return,pk_lin
end

