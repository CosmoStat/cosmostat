function mk_fisher_lens,fid,sv,cl_info=cl_info,cosmo=cosmo,cl_lens=cl_lens,cl_cov=cl_cov,dcldp_om=dcldp_om,dcldp_w0=dcldp_w0,dcldp_wa=dcldp_wa,dcldp_h=dcldp_h,dcldp_sig8=dcldp_sig8,dcldp_ob=dcldp_ob,dcldp_n=dcldp_n,dcldp_ol=dcldp_ol

; Aug 08 - modified by AA changed header and swapped Omb and sig8 order
; Aug 08 - modified by TK to allow curved model with the new fid. structure
; Aug 08 - modified by AA allow user to input derivatives (useful
;          feature for the website)
; Jul 08 - modified by AA check that sv can do lensing
; Jul 08 - modified by AA compatible with new set_fiducial & mk_cosmo
; Written by Adam Amara May 2008
; replaces an earlier routine originally written by Alex Refregier (Jan 2002)
; ***************// HELP BEGIN //**************
; PURPOSE: This routine calculates a Fisher matrix for lensing. 
;
; INPUTS: fid - Structure with fiducial parameters (this is
;                    created by the routine set_fiducial.pro)
;         sv - Structure containing survey information (this is
;              created by the routine mk_survey.pro)
;
; OPTIONAL INPUTS: cosmo - structure containing the central cosmology
;                          parameter (created by the routine
;                          mk_cosmo.pro) 
;                  cl_lens - lensing correlation functions for central
;                            cosmology (created by the routine
;                            mk_cl_tomo.pro) 
;                  cl_cov - cl covariance (created by the routine
;                           mk_cl_cov_tomo.pro)
;                  dcld? - If you have pre-calculated the derivatives
;                          these can be input to speed up the calculation 
; OUTPUTS: Fisher - Structure containing Fisher matrix and errors
; OPTIONAL OUTPUTS: Cl_info - Structure with background information: 
;                             cosmo,cl,dcldp,pname and cl_cov
;----
; Example 1: Fisher=fisher=mk_fisher_lens(fid,sv)
; Calculates the lensing fisher matrix for the fiducial cosmology in
; the structure fiducial (created by set_fiducial) and the survey in
; the structure survey (created using mk_survey). 
;
; Example 2: Fisher=mk_fisher_lens(fisher=mk_fisher_lens(fid,sv,cosmo=cosmo,cl_lens=cl,cl_cov=cl_cov)
; Similar to example 1 but here the central cosmology and lensing
; powerspectra are input which helps to (slightly) speed up the
; calculation. 
;
; Example 3: Fisher=mk_fisher_lens(fid,sv,cl_info=cl_info)
; Similar to example 1 but the structure cl_info is output. This for
; example contains the derivatives of Cl with cosmological parameters
; ----
; ***************// HELP END //**************

; -- Start: check sv can do lensing --
void=where(sv.probes eq 'lens',count)
if (count eq 0) then begin
   print,'Error:This survey cannot be used for gravitational lensing'
   return,-1
endif
; -- End: check sv can do lensing --

; -- Start: Calculate Cl and Cl_cov for central model --
if not keyword_set(cosmo) then cosmo=mk_cosmo(fid) 
if not keyword_set(cl_lens) then cl_lens=mk_cl_tomo(fid,cosmo,sv) 
if not keyword_set(cl_cov) then cl_cov=mk_cl_cov_tomo(fid,cl_lens,sv) 
; -- End: Calculate Cl and cl_cov for central model --

; -- Start: Define Fisher matrix parameters --
; construct vector of central values parameters:
p=[cosmo.const.omega_m,cosmo.const.w0,cosmo.const.wa,$
   cosmo.const.h,cosmo.const.omega_b,cosmo.const.sigma8,$
   cosmo.const.n,cosmo.const.omega_l]
pname=['!7X!6!im!n','w0','wa','h','!7X!6!ib!n','!4o!6!d8!n','n','!7X!6!il!n']
;Check that the number of parameters of the two parameter vectors match
if (n_elements(p) ne n_elements(pname)) then begin
   print,'The number of parameters in pname and p vectors do not match.'
   return,'Error: # of parameters.'
endif else n_p=n_elements(p)
;print cosmology parameters
if tag_check(keywords,'printit',/outvalue) then begin
   print,'Central values for cosmological parameters:'
   print,'pname:',pname
   print,'values:',p
endif
; -- End: Define Fisher matrix parameters --

; -- Start: Calculate derivatives --
; array that will contain dCl/dp functions
n_cl=cl_lens.n_cl
n_l=fid.calc.n_l
dcldp=dblarr(n_l,n_p,n_cl)

; vary omega_m
if not keyword_set(dcldp_om) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_om=mk_dcldp_2p(fid,sv,'omega_m') else $
         dcldp_om=mk_dcldp_4p(fid,sv,'omega_m')
endif
dcldp(*,0,*)=dcldp_om

; vary w_0
if not keyword_set(dcldp_w0) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_w0=mk_dcldp_2p(fid,sv,'w0') else $
         dcldp_w0=mk_dcldp_4p(fid,sv,'w0') 
endif
dcldp(*,1,*)=dcldp_w0

; vary w_a
if not keyword_set(dcldp_wa) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_wa=mk_dcldp_2p(fid,sv,'wa') else $
         dcldp_wa=mk_dcldp_4p(fid,sv,'wa') 
endif
dcldp(*,2,*)=dcldp_wa

; vary h
if not keyword_set(dcldp_h) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_h=mk_dcldp_2p(fid,sv,'h') else $
         dcldp_h=mk_dcldp_4p(fid,sv,'h') 
endif
dcldp(*,3,*)=dcldp_h

; vary omega_b
if not keyword_set(dcldp_ob) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_ob=mk_dcldp_2p(fid,sv,'omega_b') else $
         dcldp_ob=mk_dcldp_4p(fid,sv,'omega_b') 
endif
dcldp(*,4,*)=dcldp_ob

; vary sigma_8
if not keyword_set(dcldp_sig8) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_sig8=mk_dcldp_2p(fid,sv,'sigma8') else $
         dcldp_sig8=mk_dcldp_4p(fid,sv,'sigma8')
endif
dcldp(*,5,*)=dcldp_sig8


; vary n
if not keyword_set(dcldp_n) then begin
   if keyword_set(fid.calc.speed) then $
      dcldp_n=mk_dcldp_2p(fid,sv,'n') else $
         dcldp_n=mk_dcldp_4p(fid,sv,'n') 
endif
dcldp(*,6,*)=dcldp_n

if (fid.cosmo.curv ge 1b) then begin
; vary omega_l
   if not keyword_set(dcldp_ol) then begin
      if keyword_set(fid.calc.speed) then $
         dcldp_ol=mk_dcldp_2p(fid,sv,'omega_l') else $
            dcldp_ol=mk_dcldp_4p(fid,sv,'omega_l') 
   endif
   dcldp(*,7,*)=dcldp_ol	
endif
;-- End: Calculate derivatives --

; For Flat cosmology need to remove Omega_v entries
if (fid.cosmo.curv ne 1b) then begin
   keep = make_array(n_p,value=1)
   ind=where(pname eq '!7X!6!il!n')
   keep(ind) = -666
   ii=where(keep eq 1)
   pname=pname(ii)
   dcldp=dcldp(*,ii,*)
   n_p=n_p-1    
endif

; compute fisher matrix
if (fid.calc.verbose ge 1) then print,'Computing Fisher Matrices'
f=dblarr(n_p,n_p)          ; all power spectra with full covariance
intgd=dblarr(n_l)
f_diag=dblarr(n_p,n_p)     ; all power spectra with only diagonal errors
f_auto=dblarr(n_p,n_p)     ; auto-correlation power spectra with diag errors (old method)
for i=0,n_p-1 do begin
  for j=0,n_p-1 do begin
    for k=0,n_l-1 do $
      intgd(k)=transpose(reform(dcldp(k,i,*)))#$
                  invert(reform(cl_cov.cov_l(k,*,*)))#reform(dcldp(k,j,*))
      f(i,j)=int_tabulated(alog(cl_lens.l),cl_lens.l*intgd,/double)
  endfor
endfor

cov=la_invert(f)
sigma_fix=dblarr(n_p) & sigma_marg=dblarr(n_p)
sigma_fix_diag=dblarr(n_p) & sigma_marg_diag=dblarr(n_p)
sigma_fix_auto=dblarr(n_p) & sigma_marg_auto=dblarr(n_p)
for i=0,n_p-1 do begin
  sigma_fix(i)=sqrt(1./f(i,i)) & sigma_marg(i)=sqrt(cov(i,i))
endfor

if keyword_set(printit) then begin
  print,'Fisher matrix: F_ij:',f
  print,'sigma (fixed):',sigma_fix
  print,'cov(pi,pj):',cov
  print,'sigma (marginalised):',sigma_marg
endif

;check that the fisher and covariance matrices have positive diagonals
;and eigen values:
f_check=check_matrix(f)
cov_check=check_matrix(cov)
status={f_d:f_check.CHECK_D,f_ev:f_check.CHECK_ev,$
        cov_d:cov_check.CHECK_D,cov_ev:cov_check.CHECK_ev}

;construct a structure with extra cl information:
Cl_info={cosmo:cosmo,cl:cl_lens,dcldp:dcldp,pname:pname,cl_cov:cl_cov}

;construct the fisher matrix:
fisher={pname:pname(0:n_p-1),p:p(0:n_p-1),f:f,cov:cov,sigma_marg:sigma_marg,sigma_fix:sigma_fix,status:status}

return,fisher
end
